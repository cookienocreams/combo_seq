module combo_seq

#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
# Combo-Seq pipeline
#
# This script can be used to analyze Combo-Seq libraries using the pseudo-aligner Salmon
#
#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
using CairoMakie
using CSV
using DataFrames
using ProgressMeter
using MultivariateStats
using StatsBase
using StatsPlots
using Statistics
using MLBase
using GLM
using Measures
using GZip
using UMAP
using Clustering
using Distances
using PDFmerger
using DataStructures
using XAM
using BGZFStreams
using CodecZlib
using ArgParse

# Define a custom exception for missing references
struct MissingReferenceError <: Exception
    message::String
end

# Define a custom exception for missing files
struct MissingFileError <: Exception
    message::String
end

# Define the Config struct to store command-line arguments
struct Config
    need_reference::Bool
    mrna::Union{String, Nothing}
    fasta::Union{String, Nothing}
    transcript::Union{String, Nothing}
    genome::Union{String, Nothing}
    organism::Union{String, Nothing}
    threads::Int

    # Constructor with error checks, Union stores possibility of no argument being passed
    function Config(need_reference::Bool, mrna::Union{String, Nothing}, fasta::Union{String, Nothing}, 
                    transcript::Union{String, Nothing}, genome::Union{String, Nothing}
                    , organism::Union{String, Nothing}, threads::Int
        )
        
        # Check if user specified a Salmon reference   
        if isnothing(mrna) && !need_reference
            throw(MissingReferenceError("No Salmon reference specified"))
        end

        # Check if user specified a miRNA reference
        if fasta == "data/mirgene_all.fas" && !isfile(fasta)
            throw(MissingFileError("No miRNA reference found"))
        end
        
        new(need_reference, mrna, fasta, transcript, genome, organism, threads)
    end
end

"""
    capture_target_files(files_to_capture::AbstractString, directory::AbstractString=".")

List all files in a directory. 

Check to see if each file contains the target file's string. 

# Example
```julia
julia> capture_target_files(".txt")
3-element Vector{String}:
 "file1.txt"
 "file2.txt"
 "file3.txt"
```
"""
function capture_target_files(files_to_capture::AbstractString, directory::AbstractString=".")
    return filter(file -> occursin(files_to_capture, file), readdir(directory))
end

"""
Update progress bar on the command line.
"""
function progress_bar_update(number_of_records::Int64
                            , interval::Number
                            , description::String
                            )
    progress_bar_update = Progress(number_of_records
    , dt=interval
    , barglyphs=BarGlyphs("[=> ]")
    , barlen=100
    , desc=description
    )

    return progress_bar_update
end

"""
    get_read_q_score!(line::String, q_score_list::Vector{Number})

Calculate the quality score for a given quality string.

This function accepts ASCII-encoded quality score and produces the average q score 
using a scalar value that estimates the 'total error probability' for that record.
It means that if you want to calculate average quality score, you don't just sum
up all the phred scores and find the average, but rather consider the 
probability distribution of the scores.

The Phred score for a base `Q` is calculated as `Q = ASCII value of quality score - 33`.
The error probability `P` for that base is then calculated as `P = 10^(-Q/10)`.
Then these probabilities are summed up for all the bases to derive total error.

# Arguments

`q_score` - The quality scores of the current line.
`q_score_list` - The vector containing all calculated quality scores.

# Returns

The updated q score list the line's average quality score.

# Reference

Illumina's explanation on quality scores and their corresponding error rates:
<https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html>

# Example
```julia
julia> get_read_q_score!("FGFEEE<FC<FGGGGGGGGGGGFFGFG8<@8CFFGF8EFGGCGFGGGGGGFG", [36.2, 35.9])
3-element Vector{Float64}:
 36.2
 35.9
 35.7
```
"""
function get_read_q_score!(line::String, q_score_list::Vector{Float64})
    probability_sum = 0

    # Find average read quality by converting quality encodings to Unicode numbers
    for char in codeunits(line)
        phred = char - 33
        error_probability = 10^(-phred / 10.0)
        probability_sum += error_probability
    end

    prob_q_score = probability_sum / length(line)
    q_score = -10.0 * log10(prob_q_score)

    # Add Q score to list
    push!(q_score_list, q_score) 

    return q_score_list
end

"""
Make miRNA bowtie2 reference to align fastqs with.
"""
function make_bowtie2_reference(fasta::String, organism_name::AbstractString)
    wait(run(pipeline(
            `bowtie2-build 
            $fasta 
            $organism_name`)
            , wait = false
            )
        )
end

"""
Make Salmon reference for mRNA alignment.
"""
function make_salmon_reference(organism_name::String
                                , transcriptome_fasta::String
                                , cores::Int
                                )
    # Remove extra header lines from gencode fasta file if one is used
    for line in eachline(transcriptome_fasta)
        if occursin("ENST", line) 
            println("Starting to make salmon reference, this could take a while...")
            wait(run(pipeline(
                    `salmon 
                    index 
                    -i data/$organism_name 
                    --transcripts gentrome.fa
                    -k 19 
                    --threads $cores
                    --decoys salmon_decoys.txt 
                    --gencode`)
                    , wait = false
                )
            )
            println("Finished making salmon reference")
            break
        else
            println("Starting to make salmon reference, this could take a while...")
            wait(run(pipeline(
                    `salmon 
                    index 
                    -i data/$organism_name 
                    --transcripts gentrome.fa
                    -k 19 
                    --threads $cores
                    --decoys salmon_decoys.txt`)
                    , wait = false
                )
            )
            println("Finished making salmon reference")
            break
        end
    end
end

"""
Function asks the user if they would like to download a non-human reference
for use by Salmon during the alignment step. If yes, the files are downloaded
and the genome names are extracted to name the Salmon directory.
"""
function create_reference_file(need_reference::Bool
                                , reference_transcriptome_fasta::String
                                , reference_fasta::String
                                , organism_name::String
                                , fasta_file::String
                                , cores::Int
                                )

    if need_reference
        # Make diredctory to store downloaded files and created reference files
        !isdir("data") && mkdir("data")

        # Check to see if the url contains the secure hypertext protocol, i.e. 'https:'
        # If yes, that is replaced by the non-secure protocol, i.e. http:, to avoid certificate errors
        transcriptome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_transcriptome_fasta)
        if transcriptome_fasta_hypertext_protocol.match == "https:"
            reference_transcriptome_fasta_address = match(r"(?!.+\/\/).+", reference_transcriptome_fasta)
            reference_transcriptome_fasta = "http:" * reference_transcriptome_fasta_address.match
        end

        # Capture filenames to be used to refer to the downloaded file
        reference_transcriptome_fasta_name = last(splitpath(reference_transcriptome_fasta))
        unzipped_reference_transcriptome_fasta = reference_transcriptome_fasta_name[1:end-3]

        # Download desired transcriptome fasta file and determine the filename
        run(`wget 
            $reference_transcriptome_fasta 
            -O $reference_transcriptome_fasta_name
            --no-check-certificate
            --quiet`
        )

        genome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_fasta)
        if genome_fasta_hypertext_protocol.match == "https:"
            reference_genome_fasta_address = match(r"(?!.+\/\/).+", reference_fasta)
            reference_fasta = "http:" * reference_genome_fasta_address.match
        end

        reference_fasta_name = last(splitpath(reference_fasta))
        unzipped_reference_fasta_name = reference_fasta_name[1:end-3]

        # Make miRNA fasta for bowtie2 reference
        create_target_organism_fasta_file(organism_name, fasta_file)

        # Download desired genome fasta file
        run(`wget 
            $reference_fasta 
            -O $reference_fasta_name
            --no-check-certificate
            --quiet`
        )

        # Unzip downloaded fasta before running it through Salmon index generation
        run(`gunzip $reference_transcriptome_fasta_name`)
        run(`gunzip $reference_fasta_name`)
        
        # Create decoy sequences list needed to mask the genome and improve mRNA quantification
        # See Salmon documentation for more info https://salmon.readthedocs.io/en/latest/salmon.html
        genomic_target_names = 
        read(pipeline(
                    `cat $unzipped_reference_fasta_name`
                    ,`grep "^>"`
                    , `cut -d " " -f 1`
                    )
                    , String
        )

        # Write decoy sequences to a new file for use during Salmon index creation
        open("salmon_decoys.txt", "w") do decoys_output
            write(decoys_output, genomic_target_names)
        end
        run(`sed -i.bak -e 's/>//g' salmon_decoys.txt`)

        # Add genome fasta to transcriptome fasta for masking during Salmon index creation 
        concatenated_reference = read(pipeline(`
                                    cat $unzipped_reference_transcriptome_fasta
                                    $unzipped_reference_fasta_name
                                            `)
        )
        open("gentrome.fa", "w") do concatenated_reference_file
            write(concatenated_reference_file, concatenated_reference)
        end

        # Create reference Salmon and Bowtie2 databases for use in further analyses
        make_salmon_reference(organism_name, unzipped_reference_transcriptome_fasta, cores)
        make_bowtie2_reference(fasta_file, organism_name)

        trash_removal("gentrome.fa")
        trash_removal("salmon_decoys.txt")
        trash_removal(unzipped_reference_fasta_name)
        trash_removal(unzipped_reference_transcriptome_fasta)

    else
        !isdir("data") && mkdir("data")

        create_target_organism_fasta_file(organism_name, fasta_file)
    end

    return organism_name
end

"""
    parse_fastq_files(fastqs::Vector{String}
                        , sample_names::Vector{SubString{String}}
                        )

Loop through fastqs to calculate read counts, quality scores, and percent dimer.

Calculate the number of reads in each fastq file.

Create dictionary with each sample name and its read count.

Calculate the percent dimer present in a given fastq.

Calculate the average Q score for a given fastq. Divide sum of Unicode converted quality 
strings by the number of quality strings.

#Example
```julia
julia> parse_fastq_files(["sample1.fastq.gz","sample2.fastq.gz","sample3.fastq.gz"]
                            , ["sample1", "sample2", "sample3"])
Dict{String, Int64} with 3 entries:
  "sample1" => 6390314
  "sample2" => 5000000
  "sample3" => 7052928
Dict{String, Float64} with 3 entries:
  "sample1" => .02
  "sample2" => .007
  "sample3" => .018
Dict{String, Number} with 3 entries:
  "sample1" => 35.7
  "sample2" => 36.1
  "sample3" => 35.9
```
"""
function parse_fastq_files(fastqs::Vector{String}
                            , sample_names::Vector{SubString{String}}
                            )

    number_of_records = length(fastqs)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Parsing raw fastq files..."
                                            )

    read_count_dict = Dict{String, Int64}()
    dimer_count_dict = Dict{String, Float64}()
    q_score_dict = Dict{String, Float64}()

    for (fastq_file, sample_name) in zip(fastqs, sample_names)
        # List to store each reads quality information
        q_score_list = Vector{Float64}()
        # Store dimer information
        dimer_count = 0

        fastq = GZip.open(fastq_file)
        sample_read_count = 1
        line_tracker = 1

        for line in eachline(fastq)
            # Ignore lines without sequence or quality information
            if line_tracker % 4 != 2 && line_tracker % 4 != 0
                line_tracker += 1
            elseif line_tracker % 4 == 2
                # Count canonical 0 bp dimer
                dimer = startswith(line, "TGGAATTCTCGGGTGCC")

                if dimer
                    dimer_count += 1
                end

                sample_read_count += 1
                line_tracker += 1
            elseif line_tracker % 4 == 0
                # Calculate average quality score
                get_read_q_score!(line, q_score_list)
                line_tracker += 1
            end  
        end

        close(fastq)

        # Add read count to dictionary for use later
        read_count_dict[sample_name] = sample_read_count

        # Calculate q score for each sample
        average_q_score = round(mean(q_score_list), digits = 2)
        q_score_dict[sample_name] = average_q_score

        percent_dimer = 100 * dimer_count / sample_read_count
        dimer_count_dict[sample_name] = percent_dimer

        next!(update_progress_bar)
    end

    return read_count_dict, dimer_count_dict, q_score_dict
end

"""
    trim_adapters(fastqs::Vector{String}
                    , sample_names::Vector{SubString{String}}
                    )

Trim the 3' adapter from each read.

Combo-Seq Read 1 Setup:\n
               5' Adapter      -          Insert        - 3' Adapter      
GATCGTCGGACTGTAGAACTCTGAACNNNN - TGTCAGTTTGTCAAATACCCCA - AAAAAAAAAA 

The bases on the 3' end are also quality trimmed if their quality score is below 20. Reads 
shorter than 16 bases or that weren't trimmed are discarded.

#Example
```julia
julia> trim_adapters(["sub_sample1.fastq","sub_sample2.fastq","sub_sample3.fastq"]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1.cut.fastq"
 "sample2.cut.fastq"
 "sample3.cut.fastq"
```
"""
function trim_adapters(fastqs::Vector{String}
                        , sample_names::Vector{SubString{String}}
                        )
    number_of_records = length(fastqs)
    update_progress_bar = progress_bar_update(number_of_records, .5, "Trimming adapters...")

    for (fastq_file, sample_name) in zip(fastqs, sample_names)
        min_length = 16

        wait(run(pipeline(
                `cutadapt 
                --cores=0 
                --discard-untrimmed 
                --quality-cutoff 20 
                --cut 4 
                --adapter A"{"10"}" 
                --output $sample_name.cut.fastq 
                --minimum-length $min_length
                $fastq_file`
                , stdout="$sample_name.cutadapt_information")
                , wait=false
                )
            )

        next!(update_progress_bar)
    end

    trimmed_fastq_files::Vector{String} = capture_target_files(".cut.fastq")

    return trimmed_fastq_files
end

"""
    align_with_salmon(trimmed_fastq_files::Vector{String}
                        , salmon_reference::AbstractString
                        , sample_names::Vector{SubString{String}}
                        , need_reference::Bool
                        )
This function aligns trimmed fASTQ files to a reference genome using Salmon.

The library type (libType) is set to stranded forward (SF) since the Combo-Seq protocol is
stranded in that direction.

# Arguments

* `trimmed_fastq_files`: A vector of strings containing the paths to the trimmed fASTQ files.
* `salmon_reference`: A string containing the path to the Salmon reference genome.
* `sample_names`: A vector of strings containing the names of the samples.
* `need_reference`: A string indicating whether or not a custom Salmon index was created.

# Returns 

The function returns a vector of strings containing the paths to the Salmon output files.
"""
function align_with_salmon(trimmed_fastq_files::Vector{String}
                        , salmon_reference::AbstractString
                        , sample_names::Vector{SubString{String}}
                        , need_reference::Bool
                        )
    number_of_records::Int64 = length(trimmed_fastq_files)
    update_progress_bar = progress_bar_update(number_of_records, .5, "Aligning reads with Salmon...")

    # Iterate over the fastq files
    for (fastq_file, sample_name) in zip(trimmed_fastq_files, sample_names)

        # Check if a custom index was created
        if need_reference
            salmon_reference = string("data/", salmon_reference)
        else
            salmon_reference = salmon_reference
        end

        # Align the fastq file to the salmon reference
        wait(run(pipeline(
                        `salmon 
                        quant
                        --index $salmon_reference
                        --libType SF
                        -r $fastq_file
                        --validateMappings 
                        --fldMean 25
                        --fldSD 5
                        --writeMappings=$sample_name.sam
                        -o ./`
                        , devnull
                        )
                        , wait = false
            )
        )

        # Rename the output files
        try
            mv("quant.sf", "$sample_name.quant.sf", force=true)
            mv("logs/salmon_quant.log", "$sample_name.salmon_quant.log", force=true)
            mv("lib_format_counts.json", "$sample_name.lib_format_counts.json", force=true)
        catch
            println("Error, no Salmon output found.")
            println("Exiting program...")
            exit()
        end

        # Remove the transcript aligned sam files
        trash_removal(capture_target_files(".mRNA.sam"))

        # Update the progress bar
        next!(update_progress_bar)
    end

    # Get the names of the salmon output files
    salmon_mrna_counts = capture_target_files("quant.sf")

    return salmon_mrna_counts
end

"""
Use aligned SAM files to create a read length distribution file.
"""
function calculate_read_length_distribution(sam_files::Vector{String}
                                            , sample_names::Vector{SubString{String}}
                                            , cores::Int
                                            )
    number_of_records = length(sam_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Calculating Read Length Distribution..."
                                            )

    for (sam_file, sample_name) in zip(sam_files, sample_names)
	    #=
        If the SAM file is empty or almost empty, it will cause errors 
        downstream so the script exits.
        =#
        read_length_distibution = 
        try
            read(pipeline(
                `samtools 
                stats 
                -@ $cores
                $sam_file`
                , `grep ^RL`
                , `cut -f 2-`)
                , String
            )
        catch
            #=
            If the SAM file is empty or almost empty, it will cause errors 
            downstream so the script exits.
            =#
            println("Error, SAM file contains no read length distribution.")
            println("Exiting program...")
            exit()
        end

        next!(update_progress_bar)

        # Write read length data to output file
        output_length_file::IOStream = open("read_lengths_$sample_name.tsv", "w")
        write(output_length_file, string("Length", "\t", "Reads", "\n"))
        write(output_length_file, read_length_distibution)

        close(output_length_file)
    end

    length_files::Vector{String} = capture_target_files("read_lengths_")

    return length_files
end

"""
Create barplot of each sample's fragment lengths based on the number of reads found at each length
"""
function plot_fragment_lengths(length_files::Vector{String}
                                , sample_names::Vector{SubString{String}}
                                )
    number_of_records = length(length_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , 1
                                            , "Creating Read Length Plot..."
                                            )

    for (file, sample_name) in zip(length_files, sample_names)
        length_file::DataFrame = CSV.read(file, DataFrame)

        # Isolate fragment lengths and the number of reads at each length
        fragment_lengths::Vector{Int64} = length_file[!, :Length]
        reads_per_length::Vector{Int64} = length_file[!, :Reads]
        next!(update_progress_bar)

        plot = Figure(resolution = (1800, 1200))
        next!(update_progress_bar)

        Axis(plot[1, 1]
            , xticks = 15:5:maximum(fragment_lengths)
            , xlabel = "Fragment Lengths"
            , ylabel = "Number of Reads"
            , title = sample_name
        )
        next!(update_progress_bar)

        # Make sample barplot
        barplot!(fragment_lengths
                , reads_per_length
                , fillto = -1
                , color = fragment_lengths
                , strokecolor = :black
                , strokewidth = 1
        )

        save(string(sample_name, "_fragment_lengths.pdf"), plot)
        next!(update_progress_bar)
    end

    pdf_length_files::Vector{String} = capture_target_files("_fragment_lengths.pdf")

    merge_pdfs([pdf_length_files...], "ComboSeq_fragment_length_plots.pdf")

    trash_removal(pdf_length_files)
end

"""
    calculate_salmon_metrics(trimmed_fastq_files::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            , dimer_count_dict::Dict{String, Float64}
                            , q_score_dict::Dict{String, Number}
                            )

Align fastq with specified reference RNA types.

Divide reads aligned with total reads to calculate percent alignment. Store 
alignment information in a dictionary.

#Example
```julia
julia> calculate_salmon_metrics("sample1.cut.fastq"
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928))
Dict{String, Float64} with 5 entries:
  "miRNA" => 65.0
  "tRNA" => 4.2
  "piRNA" => 1.2
  "snoRNA" => 1.9
  "rRNA" => 3.3
```
"""
function calculate_salmon_metrics(trimmed_fastq_files::Vector{String}
                                    , read_count_dict::Dict{String, Int64}
                                    , sample_names::Vector{SubString{String}}
                                    , dimer_count_dict::Dict{String, Float64}
                                    , q_score_dict::Dict{String, Float64}
                                    )
    number_of_records = length(trimmed_fastq_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Calculating Sample Metrics..."
                                            )

    # Make output metrics file to write stats into
    metrics_file = open("ComboSeq_metrics.csv", "w")

    # Write metrics header information
    write(metrics_file, string(
        "Sample", ","
        , "Read Count", ","
        , "Aligned Reads", ","
        , "Average Q score", ","
        , "Dimer", ","
        , "<16 bp", ","
        , "miRNA", ","
        , "mRNA", ","
        , "Directionality", ","
        , "Unaligned Reads"
        , "\n"
        )
    )

    for sample_name in sample_names
        # Calculate q_score score
        average_q_score = q_score_dict[sample_name]

        # Gather the number of canonical dimer reads
        percent_dimer = dimer_count_dict[sample_name]

        # Gather metrics data from accumlulated analysis files
        cutadapt_information = read(`cat $sample_name.cutadapt_information`, String)
        salmon_library_counts_file = read(`cat $sample_name.lib_format_counts.json`, String)
        salmon_log_file = read(`cat $sample_name.salmon_quant.log`, String)
        miRNA_aligned_reads = readchomp(`samtools view -c -F 0x4 $sample_name.miRNA.sam`)

        # Uses Salmon output JSON file which reports the number of fragments that had 
        # at least one mapping compatible with a stranded library
        total_directional_RNAs_string = match(r"\"num_assigned_fragments\":.\K\d+", salmon_library_counts_file).match
        consistent_directional_RNAs_string = match(r"\"num_frags_with_concordant.+\":.\K\d+", salmon_library_counts_file).match
        
        # Uses Salmon output log file to gather the mRNA mapping rate
        mRNA_mapping_string = match(r"Mapping rate =.\K\d+\.\d{2}", salmon_log_file).match
        
        # Search cutadapt trimming information so that adapter and short fragment data can be easily captured by Regex
        short_fragments_string = match(r"Reads that were too short:.+\(\K.+(?=\%)", cutadapt_information).match
        
        # Convert gathered strings into floats for calculations and adding to the output file
        # They need to be floats in order to plot them later
        miRNA_aligned_reads = parse(Float64, miRNA_aligned_reads)
        short_fragments = abs(parse(Float64, short_fragments_string) - percent_dimer)
        mRNA_mapping = parse(Float64, mRNA_mapping_string)
        total_directional_RNAs = parse(Float64, total_directional_RNAs_string)
        consistent_directional_RNAs = parse(Float64, consistent_directional_RNAs_string)

        # Calculate strand directionality
        percent_directionality::Float64 = round(100 * consistent_directional_RNAs / total_directional_RNAs, digits = 2)

        # Calculate miRNA alignment rate
        miRNA_mapping::Float64 = round(100 * miRNA_aligned_reads / read_count_dict[sample_name], digits = 2)

        aligned_reads = round(miRNA_mapping + mRNA_mapping, digits = 2)
        unaligned_reads = round(100 - aligned_reads - percent_dimer - short_fragments, digits = 2)

        # Write metrics data to output file
        write(metrics_file, string(
            sample_name, ","
            , read_count_dict[sample_name], ","
            , aligned_reads, ","
            , average_q_score, ","
            , round(percent_dimer, digits = 2), ","
            , round(short_fragments, digits = 2), ","
            , miRNA_mapping, ","
            , mRNA_mapping, ","
            , percent_directionality, ","
            , unaligned_reads
            , "\n"
            )
        )

        next!(update_progress_bar)
    end

    close(metrics_file)

    # Returns the metrics file name
    return "ComboSeq_metrics.csv"
end

"""
Create violin plot of each sample's RNA aligment metrics
"""
function make_metrics_violin_plot(metrics_file::String)

    # Import sample metrics file
    metrics_file::DataFrame = CSV.read(metrics_file, DataFrame)

    # Create new dataframe without sample names and add sample and column names to arrays
    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    number_of_records = length(column_names)
    update_progress_bar = progress_bar_update(number_of_records, .25
                                            , "Creating Violin Plots..."
                                            )
    
    # Gather each metric's column data for plotting
    metrics_columns(a) = [number for number in metrics_file[!, a]]

    # Make empty figure and set x and y axes
    # x-axis is vector of integers from 1 to the number of metrics
    # y-axis is vector of vectors containing the data from each metric
    fig = Figure(resolution = (1800, 1200))
    x, y = 1:num_of_cols, [metrics_columns(column) for column in 1:num_of_cols]
    
    # Create violin plot with data from all samples analyzed
    for column in 1:num_of_cols
        ax = Axis(fig[1, column]
                , yticks = 0:5:100
                , ylabel = "Percent"
                , title = column_names[column]
        )
        CairoMakie.ylims!(ax, 0, 100)

        # Make violin plot with combined sample data
        CairoMakie.violin!(fig[1, column]
                            , repeat([x[column]]
                            , first(size(metrics_file))
                            )
                            , y[column]
                            , show_median=true
        )
        next!(update_progress_bar)

        # Save file once all columns have been added
        if column == num_of_cols
            save("Violin_plot_metrics.png", fig)
        end

        next!(update_progress_bar)
    end
end

"""
Create barplot of each sample's RNA aligment metrics.
"""
function plot_metrics(metrics_file::String
                        , sample_names::Vector{SubString{String}}
                        )

    metrics_file::DataFrame = CSV.read(metrics_file, DataFrame)

    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    colors = [:snow3, :honeydew2, :lightpink, :bisque2, :deepskyblue, :aquamarine3
            , :peachpuff, :paleturquoise1, :orange1, :lightsalmon2, :lemonchiffon
    ]
    
    sample_rows(a) = [values for values in metrics_file[a, :]]

    number_of_records = length(sample_names)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .25
                                            , "Creating Metrics Plots..."
                                            )
    
    for sample in eachindex(sample_names)
        sample_name = sample_names[sample]

        # Set axis values; y-axis is a vector containing each sample's calculated metrics
        x, y = 1:num_of_cols, sample_rows(sample)
        fig = Figure(resolution = (1800, 1200))

        ax = Axis(fig[1, 1]
                , xticks = (1:num_of_cols
                , column_names)
                , yticks = 0:10:100
                , ylabel = "Percent"
                , title = string(sample_name, " Metrics - ", "Combo-Seq")
        )

        # Set y-axis limits
        CairoMakie.ylims!(ax, 0, 105)

        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , strokecolor = :black
                , color = colors[1:num_of_cols]
                , bar_labels = :y
        )

        next!(update_progress_bar)

        save(string(sample_name, "_metrics.pdf"), fig)
    end

    pdf_metrics_files::Vector{String} = capture_target_files("_metrics.pdf")

    merge_pdfs([pdf_metrics_files...], "ComboSeq_metrics_plots.pdf")
    
    trash_removal(pdf_metrics_files)
end

"""
Create a file with a list of unique species names given a miRNA mature or hairpin fasta.
"""
function get_genus_names(fasta::String)
    genuses = Set{String}()

    open(fasta) do miRNA_reference
        for line in eachline(miRNA_reference)
            if startswith(line, ">")
                genus_abbreviation = line[2:4]
                push!(genuses, genus_abbreviation)
            end
        end
    end

    return genuses
end

"""
Function takes in mirBase mature miRNA fasta file with data from all available
organisms and pulls out only the miRNA data pertaining to the target organism.
"""
function create_target_organism_fasta_file(organism_name::String, fasta_file::String)

    # Skip making species miRNA reference if it exists in data/
    bowtie2_reference_files = capture_target_files(".bt2", "data")
    if isempty(bowtie2_reference_files)

        # Create file with all genuses with miRNA annotations
        genuses = get_genus_names(fasta_file)

        # Set the target or organism's genus and species for labeling the output file
        target_genus_abbreviation = ""
        output_fasta_file_name = organism_name * ".fa"
        bowtie2_reference_name = organism_name

        # Loop through file containing the genus of all the organisms with miRNA data
        # If the target organism is in that list, the genus abbreviation is stored for use below
        for genus in genuses
            if occursin(lowercase(organism_name), lowercase(genus))
                target_genus_abbreviation = genus
                break
            end
        end

        # Open the IO streams and sets the output file names for the fasta and bowtie2 index
        input_fasta_file = open(fasta_file)
        output_fasta_file = open(string("data/", organism_name, ".fa"), "w")

        #=
        Loop through miRNA fasta looking for the header and sequence information
        for the target genus. If there's a match, that information is added to a new fasta file.
        Must convert RNA sequences in miRNA fasta to DNA.
        =#
        for line in eachline(input_fasta_file)
            input_genus_abbreviation = match(r">\K\w+", line)
            if !isnothing(input_genus_abbreviation)
                if target_genus_abbreviation == input_genus_abbreviation.match
                    write(output_fasta_file, line, "\n")
                    write(output_fasta_file, string(replace(readline(input_fasta_file), "U" => "T"), "\n"))
                end
            end
        end

        close(input_fasta_file)
        close(output_fasta_file)

        #=
        Use the newly created single organism fasta to build a bowtie2 index; will be used
        for alignment with bowtie2 to determine miRNA count information.
        =#
        wait(run(pipeline(
            `bowtie2-build 
            data/$output_fasta_file_name 
            data/$bowtie2_reference_name`
            , devnull)
            , wait = false
            )
        )
        
    end
end

"""
    generate_mirna_counts(input_sam_file::AbstractString)

This function creates a DataFrame containing all aligned miRNA and the number of times 
they appear.

# Returns

A DataFrame with the following columns:

* `name`: The name of the miRNA.
* `count`: The number of times the miRNA appears in the SAM file.
"""
function generate_mirna_counts(input_sam_file::AbstractString)
    miRNA_counts = Dict{String, Int64}()

    # Open the input sam file
    sam_file = open(SAM.Reader, input_sam_file)
    record = SAM.Record()

    # Go through each sequence in the sam file
    while !eof(sam_file)
        # Empty sam record to reduce memory usage
        empty!(record)
        read!(sam_file, record)

        # Unmapped reads have no ID
        if SAM.ismapped(record)
            miRNA_name = SAM.refname(record)

            # Check if the miRNA has already been added, update count if needed
            if haskey(miRNA_counts, miRNA_name)
                miRNA_counts[miRNA_name] += 1
            else
                miRNA_counts[miRNA_name] = 1
            end
        end
    end

    # Create a DataFrame with the miRNA names and counts
    miRNA_counts_dataframe = DataFrame(name = collect(keys(miRNA_counts))
                                    , count = collect(values(miRNA_counts))
                                    )
    # Sort the DataFrame by miRNA count in descending order
    sort!(miRNA_counts_dataframe, :count, rev = true)

    return miRNA_counts_dataframe
end

"""
    mirna_discovery_calculation(trimmed_fastq_files::Vector{String}
                                , sample_names::Vector{SubString{String}}
                                , bowtie2_reference_name::AbstractString
                                , organism_name::AbstractString
                                , need_reference::AbstractString
                                )

Align trimmed fastq files to the single organism miRNA bowtie2 reference. 

The bowtie2 output SAM is input into function to generate miRNA read counts.

#Example
```julia
julia> mirna_discovery_calculation(["sample1.cut.fastq","sample2.cut.fastq","sample3.cut.fastq"]
                                , ["sample1", "sample2", "sample3"])
3-element Vector{DataFrame}:
469x2 DataFrame
 Row │ name             count
     │ String           Int64
─────┼────────────────────────
   1 │ hsa-miR-1307-3p      1
   2 │ hsa-miR-425-5p      58
   3 │ hsa-miR-488-5p       2
 [...]

 3-element Vector{String}:
 "sample1.miRNA.sam"
 "sample2.miRNA.sam"
 "sample3.miRNA.sam"
```
"""
function mirna_discovery_calculation(trimmed_fastq_files::Vector{String}
                                    , sample_names::Vector{SubString{String}}
                                    , reference::AbstractString
                                    , cores::Int
                                    )
    number_of_records = length(trimmed_fastq_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Calculating the number of miRNA present..."
                                            )

    vector_of_miRNA_counts_dfs = Vector{DataFrame}()

    for (fastq_file, sample_name) in zip(trimmed_fastq_files, sample_names)

        # Align to miRNA reference
        wait(run(pipeline(
                `bowtie2
                --norc
                --threads $cores
                -x data/$reference
                -U $fastq_file
                -S $sample_name.miRNA.sam`
                , devnull)
                , wait = false
                )
            )

        counts_df = generate_mirna_counts(string(sample_name, ".miRNA.sam"))
        push!(vector_of_miRNA_counts_dfs, counts_df)        
        next!(update_progress_bar)
    end

    sam_files::Vector{String} = capture_target_files(".miRNA.sam")

    return vector_of_miRNA_counts_dfs, sam_files
end

"""
    plot_mirna_counts(mirna_counts_dfs::Vector{DataFrame}
                        , sample_names::Vector{SubString{String}}
                        )

Plot counts of the number of unique miRNA each sample aligned to.

#Example
```julia
julia> plot_mirna_counts([sample1_counts_df,sample2_counts_df,sample3_counts_df]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1_miRNA_counts.csv"
 "sample2_miRNA_counts.csv"
 "sample3_miRNA_counts.csv"
```
"""
function plot_mirna_counts(mirna_counts_dfs::Vector{DataFrame}
                            , sample_names::Vector{SubString{String}}
                            )
    number_of_records = length(mirna_counts_dfs)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .25
                                            , "Counting the miRNA in each sample..."
                                            )

    for (miRNA_df, sample_name) in zip(mirna_counts_dfs, sample_names)
        # Calculate the total number of reads mapped to miRNA
        total_mapped_reads = sum(miRNA_df[!, :count])

        # Determine reads per million (RPM) for each sample
        RPM::DataFrame = combine(miRNA_df
        , :count => ByRow(miRNA_reads -> round(miRNA_reads / total_mapped_reads * 10^6, digits = 2)
        ))

        # Remove "(unk)" from all miRNA names
        clean_miRNA_names::DataFrame = combine(miRNA_df
                                                , :name => ByRow(name -> replace(name, "(unk)" => ""))
        )

        # Reconstitute dataframe with new miRNA names, miRNA counts, and RPM column
        miRNA_df::DataFrame = hcat(clean_miRNA_names, select(miRNA_df, :count), RPM)

        # Rename column with miRNA names and RPM data
        rename!(miRNA_df, :count_function => :RPM, :name_function => :name)

        # Check column containing miRNA counts and adds miRNA above a set threshold
        threshold_of_one = length(filter(>=(1), miRNA_df[!, :count]))
        threshold_of_three = length(filter(>=(3), miRNA_df[!, :count]))
        threshold_of_five = length(filter(>=(5), miRNA_df[!, :count]))
        threshold_of_ten = length(filter(>=(10), miRNA_df[!, :count]))

        # Barplot colors
        colors = [:grey88, :skyblue2, :peachpuff, :lightsalmon]

        # Set x and y-axis values; y-axis is a vector containing each sample's counted miRNAs
        x, y = 1:4, [threshold_of_one, threshold_of_three, threshold_of_five, threshold_of_ten]

        fig = Figure(resolution = (1800, 1200))
        increments = round(threshold_of_one / 10, sigdigits=2)
        ax = Axis(fig[1, 1]
                , xticks = (1:4
                , ["Threshold 1", "Threshold 3", "Threshold 5", "Threshold 10"])
                , yticks = 0:increments:threshold_of_one
                , ylabel = "Number of miRNA"
                , title = string(sample_name, " miRNA Counts")
        )
        CairoMakie.ylims!(ax, 0, threshold_of_one * 1.05)

        # Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
        )

        save(string(sample_name, "_miRNA_counts.pdf"), fig)

        # Write output file containing miRNA names, read count, and RPM data
        CSV.write(string(sample_name, "_miRNA_counts.csv"), miRNA_df)
        
        next!(update_progress_bar)
    end

    pdf_miRNA_files::Vector{String} = capture_target_files("_miRNA_counts.pdf")

    # Merge output PDFs into one file
    merge_pdfs([pdf_miRNA_files...], "Comboseq_miRNA_plots.pdf")
    
    # Remove individual sample plots
    trash_removal(pdf_miRNA_files)

    miRNA_counts_files::Vector{String} = capture_target_files("_miRNA_counts.csv")

    return miRNA_counts_files
end

function plot_mrna_counts(salmon_mrna_counts::Vector{String}
                            , sample_names::Vector{SubString{String}}
                            )
    number_of_records = length(salmon_mrna_counts)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .25
                                            , "Counting the mRNA in each sample..."
                                            )

    for (file, sample_name) in zip(salmon_mrna_counts, sample_names)
        # Import Salmon output file containing mRNA information
        mRNA_file::DataFrame = CSV.read(file, DataFrame)

        # Checks column containing mRNA counts and adds mRNA above a set threshold
        threshold_of_one = sum(mRNA_file[!,:TPM] .>= 1)
        threshold_of_ten = sum(mRNA_file[!,:TPM] .>= 10)
        threshold_of_one_hundred = sum(mRNA_file[!,:TPM] .>= 100)
        threshold_of_five_hundred = sum(mRNA_file[!,:TPM] .>= 500)
        threshold_of_one_thousand = sum(mRNA_file[!,:TPM] .>= 1000)

        # Captures the name and TPM count of each mRNA with a TPM of at least 1
        mRNA_at_threshold_of_one = mRNA_file[mRNA_file[!,:TPM] .>= 1, [:Name, :TPM]]

        # Sort output dataframe with mRNA TPM counts by the most highly expressed mRNAs
        sorted_mRNA_at_threshold_of_one = sort(mRNA_at_threshold_of_one[!, [:Name, :TPM]]
                                                , :TPM
                                                , rev = true
        )

        # Barplot colors
        colors = [:burlywood, :lightsteelblue, :mistyrose2, :pink4, :grey80]

        # Set x and y-axis values; y-axis is a vector containing each sample's counted mRNAs 
        x = 1:5
        y = [threshold_of_one
            , threshold_of_ten
            , threshold_of_one_hundred
            , threshold_of_five_hundred
            , threshold_of_one_thousand
        ]

        fig = Figure(resolution = (1800, 1200))
        increments = round(threshold_of_one / 10, sigdigits=2)
        ax = Axis(fig[1, 1]
                , xticks = (1:5
                , ["Threshold 1", "Threshold 10", "Threshold 100", "Threshold 500", "Threshold 1000"])
                , yticks = 0:increments:threshold_of_one
                , ylabel = "Number of mRNA"
                , title = string(sample_name, " mRNA Counts")
        )
        CairoMakie.ylims!(ax, 0, threshold_of_one * 1.05)

        # Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
        )

        save(sample_name * "_mRNA_counts.pdf", fig)

        # Write output file containing mRNA names and TPM data
        CSV.write(sample_name * "_mRNA_counts.csv", sorted_mRNA_at_threshold_of_one)
        
        next!(update_progress_bar)
    end

    pdf_mRNA_files::Vector{String} = capture_target_files("_mRNA_counts.pdf")

    # Merge output PDFs into one file
    merge_pdfs([pdf_mRNA_files...], "Comboseq_mRNA_plots.pdf")
    
    # Remove individual sample plots
    trash_removal(pdf_mRNA_files)

    mRNA_counts_files::Vector{String} = capture_target_files("_mRNA_counts.csv")

    return mRNA_counts_files
end

"""
    find_common_mirnas(mirna_counts_files::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            )

Return the RNAs present in all samples, together with their counts and reads per million reads (RPM)
or transcripts per million (TPM).

# Arguments
- `[mi/m]rna_counts_files`: A vector of files containing sample RNA counts.
- `read_count_dict`: A dictionary where keys are sample names and values are the total number of reads 
for the corresponding sample.
- `sample_names`: A vector of sample names.

# Returns
- A Tuple of three items:
    1. A list of RNAs that are common to all samples.
    2. A list of dictionaries where each dictionary contains the RNA counts for one sample.
    3. A list of dictionaries where each dictionary contains the RPM/TPM values for the miRNAs in one sample.

# Example
```julia
julia> find_common_mirnas(["sample1_miRNA_counts.csv", "sample2_miRNA_counts.csv"]...,
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928)
                        , ["sample1", "sample2", "sample3"])
```
"""
function find_common_rnas(mirna_counts_files::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            )
    # Vector of dictionaries containing the miRNA counts from each sample
    miRNA_info = Vector{Dict{String, Int64}}()
    RPM_info = Vector{Dict{String, Number}}()
    full_miRNA_names_list = Vector{String}()

    for (file, sample_name) in zip(mirna_counts_files, sample_names)
        miRNA_file = open(file, "r")

        # Dictionaries to hold each miRNA and its associated read count
        sample_miRNA_counts_dictionary = Dict{String, Int64}()
        sample_miRNA_RPM_dictionary = Dict{String, Number}()

        # Skip header line
        readline(miRNA_file)

        for line in eachline(miRNA_file)
            split_line = split(line, ",")
            miRNA_count = parse(Int64, split_line[2])
            miRNA_name = first(split_line)
            RPM = round(miRNA_count / (read_count_dict[sample_name] / 10^6), digits = 4)

            sample_miRNA_counts_dictionary[miRNA_name] = miRNA_count
            sample_miRNA_RPM_dictionary[miRNA_name] = RPM
        end

        push!(miRNA_info, sample_miRNA_counts_dictionary)
        push!(RPM_info, sample_miRNA_RPM_dictionary)
        
        append!(full_miRNA_names_list, keys(sample_miRNA_counts_dictionary))

        close(miRNA_file)
    end

    # Count number each miRNA's occurances
    miRNA_names_dict = countmap(full_miRNA_names_list)

    # Find miRNA present in all samples, .i.e. have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(sample_names), miRNA_names_dict) |> keys

    return miRNAs_in_common, miRNA_info, RPM_info
end

function find_common_rnas(mrna_counts_files::Vector{String}
                            , sample_names::Vector{SubString{String}}
                            )
    # Vector of dictionaries containing the mRNA counts from each sample
    mRNA_info = Vector{Dict{String, Float64}}()
    TPM_info = Vector{Dict{String, Number}}()
    full_mRNA_names_list = Vector{String}()

    for file in mrna_counts_files
        mRNA_file::IOStream = open(file, "r")

        # Dictionaries to hold each mRNA and its associated read count
        sample_mRNA_counts_dictionary = Dict{String, Float64}()
        sample_mRNA_TPM_dictionary = Dict{String, Number}()

        # Skip header line
        readline(mRNA_file)

        for line in eachline(mRNA_file)
            split_line = split(line, ",")
            if length(split_line) == 2
                mRNA_count = parse(Float64, split_line[2])
                mRNA_name = first(split_line)
                TPM = round(mRNA_count, digits = 4)

                sample_mRNA_counts_dictionary[mRNA_name] = mRNA_count
                sample_mRNA_TPM_dictionary[mRNA_name] = TPM
            end
        end

        push!(mRNA_info, sample_mRNA_counts_dictionary)
        push!(TPM_info, sample_mRNA_TPM_dictionary)
        
        append!(full_mRNA_names_list, keys(sample_mRNA_counts_dictionary))

        close(mRNA_file)
    end

    # Count number each mRNA's occurances
    miRNA_names_dict = countmap(full_mRNA_names_list)

    # Find miRNA present in all samples, .i.e. have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(sample_names), miRNA_names_dict) |> keys

    return miRNAs_in_common, mRNA_info, TPM_info
end

"""
Write RNAs all samples have in common to output file.
"""
function write_common_rna_file(mirna_names::Base.KeySet{String, Dict{String, Int64}}
                                , mirna_info::Vector{Dict{String, Int64}}
                                , rpm_info::Vector{Dict{String, Number}}
                                , sample_names::Vector{SubString{String}}
                                ; rna_type="miRNA"
                                )

    for (index, sample) in enumerate(sample_names)
        common_miRNA_file::IOStream = open(string(sample, "_common_miRNAs.tsv"), "w")
        common_miRNA_file_RPM::IOStream = open(string(sample, "_common_miRNAs_RPM.tsv"), "w")

        if index == 1
            write(common_miRNA_file, string(rna_type, "\t", sample, "\n"))
            write(common_miRNA_file_RPM, string(rna_type, "\t", sample, "\n"))
            for miRNA in mirna_names
                write(common_miRNA_file, string(miRNA, "\t", mirna_info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(miRNA, "\t", rpm_info[index][miRNA], "\n"))
            end
        else
            write(common_miRNA_file, string(sample, "\n"))
            write(common_miRNA_file_RPM, string(sample, "\n"))
            for miRNA in mirna_names
                write(common_miRNA_file, string(mirna_info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(rpm_info[index][miRNA], "\n"))
            end
        end

        close(common_miRNA_file)
        close(common_miRNA_file_RPM)
    end

    common_files = capture_target_files("_common_miRNAs.tsv")
    common_RPM_files = capture_target_files("_common_miRNAs_RPM.tsv")

    # Combine all the separate common miRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_miRNAs.tsv"))
    run(pipeline(`paste $common_RPM_files`, stdout = "Common_RPM_miRNAs.tsv"))

    # Run principal component analysis and UMAP on common miRNA
    plot_clustering("Common_RPM_miRNAs.tsv", rna_type)

    trash_removal(common_files)
    trash_removal(common_RPM_files)
end

function write_common_rna_file(mrna_names::Base.KeySet{String, Dict{String, Int64}}
                                , mrna_info::Vector{Dict{String, Float64}}
                                , tpm_info::Vector{Dict{String, Number}}
                                , sample_names::Vector{SubString{String}}
                                ; rna_type="mRNA"
                                )
    
    for (index, sample) in enumerate(sample_names)
        common_mRNA_file::IOStream = open(string(sample, "_common_mRNAs.tsv"), "w")
        common_mRNA_file_TPM::IOStream = open(string(sample, "_common_mRNAs_TPM.tsv"), "w")

        if index == 1
            write(common_mRNA_file, string(rna_type, "\t", sample, "\n"))
            write(common_mRNA_file_TPM, string(rna_type, "\t", sample, "\n"))
            for mRNA in mrna_names
                write(common_mRNA_file, string(mRNA, "\t", mrna_info[index][mRNA], "\n"))
                write(common_mRNA_file_TPM, string(mRNA, "\t", tpm_info[index][mRNA], "\n"))
            end
        else
            write(common_mRNA_file, string(sample, "\n"))
            write(common_mRNA_file_TPM, string(sample, "\n"))
            for mRNA in mrna_names
                write(common_mRNA_file, string(mrna_info[index][mRNA], "\n"))
                write(common_mRNA_file_TPM, string(tpm_info[index][mRNA], "\n"))
            end
        end

        close(common_mRNA_file)
        close(common_mRNA_file_TPM)
    end

    common_files = capture_target_files("_common_mRNAs.tsv")
    common_TPM_files = capture_target_files("_common_mRNAs_TPM.tsv")

    # Combine all the separate common mRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_mRNAs.tsv"))
    run(pipeline(`paste $common_TPM_files`, stdout = "Common_TPM_mRNAs.tsv"))

    # Run principal component analysis and UMAP on common mRNA
    plot_clustering("Common_TPM_mRNAs.tsv", rna_type)

    trash_removal(common_files)
    trash_removal(common_TPM_files)
end

"""
    plot_clustering(Common_miRNA_File::String)

Generates Principal Component Analysis (PCA) and Uniform Manifold Approximation and 
Projection (UMAP) plots for libraries given a common RNA file.

This function performs the following steps:
1. Reads the common RNA file and creates a DataFrame.
2. Extracts sample names and [mi/m]RNA names.
3. Ensures that there are at least two samples and more than one miRNA present in the dataset.
4. Transforms the data using PCA and prepares it for plotting.
5. If there are at least two principal components, plots the PCA and saves it as "common_[mi/m]RNA_PCA.png".
6. If there are at least three principal components, performs UMAP dimensionality reduction and clustering using K-medoids.
7. Plots the UMAP and saves it as "common_[mi/m]RNA_UMAP.png".
8. Writes the PCA and UMAP information to separate CSV files for cluster tracking.

# Arguments
- `Common_[mi/m]RNA_File::String`: Path to the common RNA RPM/TPM file.

# Outputs
- Creates and saves PCA and UMAP plots as "common_[mi/m]RNA_PCA.png" and "common_[mi/m]RNA_UMAP.png".
- Writes PCA and UMAP information to "PCA_information.csv" and "UMAP_information.csv".
"""
function plot_clustering(common_rna_file::String, rna_type::String)
    common_RNA_counts = DataFrame(CSV.File(common_rna_file))
    sample_names = DataFrames.names(common_RNA_counts, Not([Symbol("$rna_type")]))
    RNA_names = common_RNA_counts[!, Symbol("$rna_type")]

    # There must be at least two samples in order to perform the principal component analysis
    # Must also have at least one RNA in common
    if last(size(common_RNA_counts)) > 2 && first(size(common_RNA_counts)) > 1
        common_RNA_counts_matrix = Matrix{Float64}(select(common_RNA_counts, Not([Symbol("$rna_type")])))

        # Make PCA covariance matrix for PCA plot
        PCA_matrix = fit(PCA, common_RNA_counts_matrix; maxoutdim = 20)
        transformed_counts = predict(PCA_matrix, common_RNA_counts_matrix)

        # Make PCA covariance matrix for UMAP plot
        UMAP_PCA_matrix = fit(PCA, common_RNA_counts_matrix'; maxoutdim = 20)
        transposed_transformed_counts = predict(UMAP_PCA_matrix, common_RNA_counts_matrix')

        # Plot PCA only if there are at least two principal components
        if first(size(transformed_counts)) >= 2
            sample_name_df = DataFrame("samples" => sample_names)
            RNA_names_df = DataFrame("$rna_type" => RNA_names)

            # Write PCA information to file for cluster tracking
            PCA_values = hcat(sample_name_df, DataFrame(transformed_counts[1:2, :]', ["PC1", "PC2"]))
            CSV.write("PCA_information.csv", PCA_values)
            
            # Plot PCA
            Plots.scatter(size = (1200, 800), dpi = 300, titlefont = (16, "Computer Modern")
                            , xlabel = "PC1", ylabel = "PC2", title = "Common $rna_type: PCA"
                            , transformed_counts[1, :], transformed_counts[2, :]
                            , left_margin = 23mm, right_margin = 8mm, bottom_margin = 8mm
                            , leg = false
            )
           savefig(string("common_", rna_type, "_PCA.png"))
        end
        
        # Plot UMAP only if there are at least three principal components
        if first(size(transposed_transformed_counts)) >= 3
            sample_name_df = DataFrame("samples" => sample_names)
            RNA_names_df = DataFrame("$rna_type" => RNA_names)

            # Create low dimensional embedding for UMAP
            near_neighbors = first(size(common_RNA_counts_matrix)) - 1
            if near_neighbors < 15
                embedding = umap(transposed_transformed_counts
                                , 2
                                ; n_neighbors = near_neighbors
                                , min_dist = 0.1
            )
            else
                embedding = umap(transposed_transformed_counts
                                , 2
                                ; n_neighbors = 15
                                , min_dist = 0.1
            )
            end

            # Write UMAP information to file for cluster tracking
            UMAP_values = hcat(RNA_names_df, DataFrame(embedding', ["UMAP1", "UMAP2"]))
            CSV.write("UMAP_information.csv", UMAP_values)

            #=
            Plot UMAP of common RNAs. Color UMAP clusters based on K-medoids clustering; 
            choose number of clusters based on the number of principal components.
            =#
            distance_matrix = pairwise(Euclidean(), embedding, embedding)
            kmedoids_cluster_colors = kmedoids(distance_matrix, first(size(transposed_transformed_counts)))
            Plots.scatter(size = (1200, 800), embedding[1, :], embedding[2, :]
                            , title="Common $rna_type: UMAP", left_margin = 13mm
                            , bottom_margin = 10mm, dpi = 300
                            , marker_z = kmedoids_cluster_colors.assignments
                            , color = :lighttest, xlabel = "UMAP1", ylabel = "UMAP2"
                            , titlefont = (16, "Computer Modern"), leg = false
                            , markersize = 9, markerstrokewidth = 0.1
            )
            savefig(string("common_", rna_type, "_UMAP.png"))
        end
    end
end

"""
Spot remove unecessary intermediate files.
"""
function trash_removal(Files_to_Delete::Vector{String})
    for file in Files_to_Delete
        rm(file)
    end
end

function trash_removal(File_to_Delete::String)
    rm(File_to_Delete)
end

"""
Remove all intermediate files.
"""
function remove_intermediate_files()
    Files_to_Delete = Set(vcat(
        capture_target_files("deduplication_statistics")
        ,capture_target_files(".bt2")
        ,capture_target_files(".miRNA.")
        ,capture_target_files("v22_cluster")
        ,capture_target_files(".sam")
        ,capture_target_files(".bam")
        ,capture_target_files("quant")
        ,capture_target_files("cut")
        )
    )

    for file in Files_to_Delete
        rm(file)
    end

    json_files = capture_target_files(".json")
    trash_removal(json_files)
end

function parse_commandline()
    arguments = ArgParseSettings(prog="Combo-Seq Analysis"
                                , description = "Basic analysis of NEXTFLEX Combo-Seq libraries."
    )

    @add_arg_table! arguments begin
        "--need-reference", "-r"
            help = "Flag for specifying if a reference needs to be downloaded or not."
            action = :store_true
        "--mrna", "-M"
            help = "The full path to the Salmon mRNA reference the samples will be aligned to, \
            e.g., /home/user/hg38_mRNA."
            arg_type = String
        "--fasta", "-f"
            help = "The full mature miRNA fasta file."
            arg_type = String
            default = "data/mirgene_all.fas"
        "--transcript", "-t"
            help = "Input website url containing the desired transcriptome reference fasta \
            if one needs to be downloaded and created. For human gencode version 44 that would be \
            'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz'.
            Only needed if '--need-reference' flag is used."
            arg_type = String
        "--genome", "-g"
            help = "Input website url containing the desired genomic reference fasta \
            if one needs to be downloaded and created. For human gencode version 44 that would be \
            'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz'.
            Only needed if '--need-reference' flag is used."
            arg_type = String
		"--organism", "-O"
            help = "Abbreviated name of organism, e.g., 'hsa' for human or 'mmu' for mouse. \
            This should match the standard three letter abbreviation found in miRNA databases \
            such as miRBase and MirGeneDB."
            arg_type = String
            default = "hsa"
        "--threads", "-p"
            help = "The number of processors to use for alignment."
            arg_type = Int
            default = 12
    end

    args = parse_args(arguments)

    # Create a Config instance from parsed arguments
    config = Config(
        get(args, "need_reference", false),
        get(args, "mrna", nothing),
        get(args, "fasta", "data/mirgene_all.fas"),
        get(args, "transcript", nothing),
        get(args, "genome", nothing),
        get(args, "organism", "hsa"),
        get(args, "threads", 12)
    )

    return config
end

function julia_main()::Cint

    # Parse command line arguments from Config struct
    config = parse_commandline()

    # Say hello Issac!
    run(`echo " "`)
    run(`echo -e "\e[0;36m&&&&&&&&&...&&&&......&&&..........&"`)
    run(`echo -e   "&&&&&&&&&...&&&&&.....&&&........&&&"`)
    run(`echo -e   "&&&...&&&...&&&.&&....&&&.......&&.&&"`)
    run(`echo -e   "&&&&&&&&&...&&&..&&...&&&......&&...&&"`)
    run(`echo -e   "&&&&&&&&&...&&&...&&..&&&.....&&&&&&&&&"`)
    run(`echo -e   "&&&...&&....&&&....&&.&&&....&&.......&&"`)
    run(`echo -e   "&&&...&&&...&&&.....&&&&&...&&.........&&"`)
    run(`echo -e   "&&&...&&&&..&&&......&&&&..&&...........&&"`)
    run(`echo -e ".................. ###"`)
    run(`echo -e ".................####*"`)
    run(`echo -e "...............*######"`)
    run(`echo -e ".............*###\e[0;32m©\e[0m\e[0;36m####"`)
    run(`echo -e "...........*#########"`)
    run(`echo -e "..........*#########."`)
    run(`echo -e ".........*##########."`)
    run(`echo -e ".........*########*###*"`)
    run(`echo -e "........*##########*####"`)
    run(`echo -e ".......*##########*..*###"`)
    run(`echo -e ".....*###########.....*"`)
    run(`echo -e "....#############"`)
    run(`echo -e "...*##*##########"`)
    run(`echo -e "...*.....#########"`)
    run(`echo -e "..........########"`)
    run(`echo -e "...........*#######"`)
    run(`echo -e "............*######*"`)
    run(`echo -e "..............*#####*"`)
    run(`echo -e "................*#####"`)
    run(`echo -e "..................*###*"`)
    run(`echo -e "....................*###"`)
    run(`echo -e ".....................*###."`)
    run(`echo -e ".....................######."`)
    run(`echo -e "..................###########"`)
    run(`echo -e ".................####*..*#####\e[0m'"`)
    run(`echo -e "\t     \e[0;36m\e[1;4mP  i  p  e  l  i  n  e\e[0m"`)
    run(`echo " "`)
    run(`echo -e "\t\e[0;31mBut\e[0m \e[1;31myou\e[0m \e[1;33mcan\e[0m \e[1;32mcall\e[0m \e[0;36mme\e[0m \e[1;35mIsaac...\e[0m "`)
    run(`echo " "`)

    need_reference = config.need_reference
    fastqs = capture_target_files("_R1_001.fastq.gz")
    sample_names = map(sample -> first(split(sample, "_")), fastqs)
    organism_name = create_reference_file(need_reference
                                            , !isnothing(config.transcript) ? config.transcript : "nothing"
                                            , !isnothing(config.genome) ? config.genome : "nothing"
                                            , config.organism
                                            , config.fasta
                                            , config.threads
                                            )
    read_count_dict, dimer_count_dict, q_score_dict = parse_fastq_files(fastqs, sample_names)
    trimmed_fastq_files = trim_adapters(fastqs, sample_names)
    mirna_counts_dfs, sam_files = mirna_discovery_calculation(trimmed_fastq_files
                                                            , sample_names
                                                            , organism_name
                                                            , config.threads
                                                            )
    mirna_counts_files = plot_mirna_counts(mirna_counts_dfs, sample_names)
    salmon_mrna_counts = align_with_salmon(trimmed_fastq_files
                                        , !isnothing(config.mrna) ? config.mrna : organism_name
                                        , sample_names
                                        , need_reference
                                        )
    mrna_counts_files = plot_mrna_counts(salmon_mrna_counts, sample_names)
    length_files = calculate_read_length_distribution(sam_files, sample_names, config.threads)
    plot_fragment_lengths(length_files, sample_names)
    metrics_file = calculate_salmon_metrics(trimmed_fastq_files
                                        , read_count_dict
                                        , sample_names
                                        , dimer_count_dict
                                        , q_score_dict
                                        )
    plot_metrics(metrics_file, sample_names)
    make_metrics_violin_plot(metrics_file)
    Full_miRNA_Names_List, mirna_info, rpm_info = find_common_rnas(mirna_counts_files
                                                                , read_count_dict
                                                                , sample_names
                                                                )
    write_common_rna_file(Full_miRNA_Names_List, mirna_info, rpm_info, sample_names, rna_type="miRNA")
    full_mrna_names_list, mrna_info, tpm_info = find_common_rnas(mrna_counts_files
                                                                , sample_names
                                                                )
    write_common_rna_file(full_mrna_names_list, mrna_info, tpm_info, sample_names, rna_type="mRNA")
    remove_intermediate_files()
    println(" ")
    println("Analysis Finished")
    println(" ")

    return 0    
end

end
