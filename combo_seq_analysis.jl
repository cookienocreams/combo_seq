module combo_seq_analysis

#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
# Combo-Seq pipeline
#
# This script can be used to analyze Combo-Seq libraries using the pseudo-aligner Salmon
# Script inputs: url for desired reference transcriptome and reference genome fasta files
# Script inputs: Name of the species to be analyzed
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

"""
    CaptureTargetFiles(Files_To_Capture::String)

List all files in the current directory. 

Check to see if each file contains the target file's string. 

#Example
```julia
julia> CaptureTargetFiles(".txt")
3-element Vector{String}:
 "file1.txt"
 "file2.txt"
 "file3.txt"
```
"""
function CaptureTargetFiles(Files_To_Capture::String)
    return [file for file in readdir() if occursin(Files_To_Capture, file)]
end

"""
Update progress bar on the command line.
"""
function ProgressBarUpdate(Number_of_Records::Int64
                            , Interval::Number
                            , Description::String
                            )
    progress_bar_update = Progress(Number_of_Records
                                    , dt=Interval
                                    , barglyphs=BarGlyphs("[=> ]")
                                    , barlen=100
                                    , desc=Description
                                    )

    return progress_bar_update
end

"""
Make transcripts bowtie2 reference to align fastqs with.
"""
function MakeBowtie2Reference(Transcriptome_Fasta::String, Salmon_Reference_Name::AbstractString)
    wait(run(pipeline(
            `bowtie2-build 
            $Transcriptome_Fasta 
            $Salmon_Reference_Name`)
            , wait = false
            )
        )
end

"""
Make Salmon reference
"""
function MakeSalmonReference(Salmon_Reference_Name::String, Transcriptome_Fasta::String)
    #Remove extra header lines from gencode fasta file if one is used
    for line in eachline(Transcriptome_Fasta)
        if occursin("ENST", line) 
            wait(run(pipeline(
                    `salmon 
                    index 
                    -i $Salmon_Reference_Name 
                    --transcripts gentrome.fa
                    -k 19 
                    --threads 12
                    --decoys salmon_decoys.txt 
                    --gencode`)
                    , wait = false
                    )
                )
            break
        else 
            wait(run(pipeline(
                    `salmon 
                    index 
                    -i $Salmon_Reference_Name 
                    --transcripts gentrome.fa
                    -k 19 
                    --threads 12
                    --decoys salmon_decoys.txt`)
                    , wait = false
                    )
                )
            break
        end
    end
end

"""
Function asks the user if they would like to download a non-human reference
for use by Salmon during the alignment step. If yes, the files are downloaded
and the genome names are extracted to name the Salmon directory.
"""
function CreateReferenceFile(Need_Reference::String)

    if Need_Reference === "y"
        #Make diredctory to store downloaded files and created reference files
        !isdir("data") && mkdir("data")

        #Download files and make reference files inside data folder
        cd("data")

        #Need to gather organism name to make species specific RNA annotations for use in alignment later 
        print("What is the scientific name of the organism you would like to use (e.g. mus musculus)?: ")
        organism_name = readline()

        print("Input website url containing the desired transcriptome reference fasta: ")
        reference_transcriptome_fasta = readline()

        #Check to see if the url contains the secure hypertext protocol, i.e. 'https:'
        #If yes, that is replaced by the non-secure protocol, i.e. http:, to avoid certificate errors
        transcriptome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_transcriptome_fasta)
        if transcriptome_fasta_hypertext_protocol.match == "https:"
            reference_transcriptome_fasta_address = match(r"(?!.+\/\/).+", reference_transcriptome_fasta)
            reference_transcriptome_fasta = "http:" * reference_transcriptome_fasta_address.match
        end

        #Capture filenames to be used to refer to the downloaded file
        reference_transcriptome_fasta_name = last(splitpath(reference_transcriptome_fasta))
        unzipped_reference_transcriptome_fasta = replace(reference_transcriptome_fasta_name, ".gz" => "")

        #Download desired transcriptome fasta file and determines the filename
        run(`wget 
            $reference_transcriptome_fasta 
            -O $reference_transcriptome_fasta_name
            --no-check-certificate
            --quiet`
            )

        print("Input website url containing the desired genomic reference fasta: ")
        reference_fasta = readline()

        genome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_fasta)
        if genome_fasta_hypertext_protocol.match == "https:"
            reference_genome_fasta_address = match(r"(?!.+\/\/).+", reference_fasta)
            reference_fasta = "http:" * reference_genome_fasta_address.match
        end

        reference_fasta_name = last(splitpath(reference_fasta))
        unzipped_reference_fasta_name = replace(reference_fasta_name, ".gz" => "")
        salmon_reference_name = first(split(reference_fasta_name, "."))

        #Make miRNA fasta for bowtie2 reference
        CreateTargetOrganismFastaFile(organism_name)

        #Download desired genome fasta file
        run(`wget 
            $reference_fasta 
            -O $reference_fasta_name
            --no-check-certificate
            --quiet`
            )

        #Unzip downloaded fasta before running it through Salmon index generation
        run(`gunzip $reference_transcriptome_fasta_name`)
        run(`gunzip $reference_fasta_name`)
        
        #Create decoy sequences list needed to mask the genome and improve mRNA quantification
        #See Salmon documentation for more info https://salmon.readthedocs.io/en/latest/salmon.html
        genomic_target_names = 
        read(pipeline(
                    `cat $unzipped_reference_fasta_name`
                    ,`grep "^>"`
                    , `cut -d " " -f 1`
                    )
                    , String
            )

        #Write decoy sequences to a new file for use during Salmon index creation
        open("salmon_decoys.txt", "w") do decoys_output
            write(decoys_output, genomic_target_names)
        end
        run(`sed -i.bak -e 's/>//g' salmon_decoys.txt`)

        #Add genome fasta to transcriptome fasta for masking during Salmon index creation 
        concatenated_reference = read(pipeline(`
                                    cat $unzipped_reference_transcriptome_fasta
                                    $unzipped_reference_fasta_name
                                            `)
                                    )
        open("gentrome.fa", "w") do concatenated_reference_file
            write(concatenated_reference_file, concatenated_reference)
        end

        #Create reference Salmon and Bowtie2 databases for use in further analyses
        MakeSalmonReference(salmon_reference_name, unzipped_reference_transcriptome_fasta)
        MakeBowtie2Reference(unzipped_reference_transcriptome_fasta, salmon_reference_name)

        cd("../")
    else
        print("Would you like to use a human or mouse reference (human/mouse)?: ")
        reference_name = readline()

        !isdir("data") && mkdir("data")
        if reference_name == "mouse"
            salmon_reference_name = "/data/analysis_files/GRCm39"
            organism_name = "Mus musculus"
            cd("data")
            CreateTargetOrganismFastaFile(organism_name)
            cd("../")
        else
            #If a custom database isn't used; script defaults to using a human reference
            salmon_reference_name = "/usr/local/scripts/RNASeq_QC/gencode.v41.2"
            organism_name = "Homo Sapiens"
            cd("data")
            CreateTargetOrganismFastaFile(organism_name)
            cd("../")
        end
    end

    return salmon_reference_name, organism_name
end

"""
    GetReadQScore!(Line::String, Q_Score_List::Vector{Number})

Calculate the quality score for a given quality string.

Convert quality string to Phred33 score and update q score vector.

#Example
```julia
julia> GetReadQScore!("FGFEEE<FC<FGGGGGGGGGGGFFGFG8<@8CFFGF8EFGGCGFGGGGGGFG", [36.2, 35.9])
3-element Vector{Float64}:
 36.2
 35.9
 35.7
```
"""
function GetReadQScore!(Line::String, Q_Score_List::Vector{Number})

    #Find read quality by converting quality encodings to Unicode numbers
    q_score = sum(codeunits(Line)) / length(Line)

    #Add Q score to list and convert to a Phred33 score
    push!(Q_Score_List, q_score - 33) 

    return Q_Score_List
end

"""
    ParseFastqFiles(Fastqs::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )

Loop through fastqs to calculate read counts, quality scores, and percent dimer.

Calculate the number of reads in each fastq file.

Create dictionary with each sample name and its read count.

Calculate the percent dimer present in a given fastq.

Calculate the average Q score for a given fastq. Divide sum of Unicode converted quality 
strings by the number of quality strings.

#Example
```julia
julia> ParseFastqFiles(["sample1.fastq.gz","sample2.fastq.gz","sample3.fastq.gz"]
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
function ParseFastqFiles(Fastqs::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )

    number_of_records = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Parsing raw fastq files..."
                                            )

    read_count_dict = Dict{String, Int64}()
    dimer_count_dict = Dict{String, Float64}()
    q_score_dict = Dict{String, Number}()

    for (fastq_file, sample_name) in zip(Fastqs, Sample_Names)
        #List to store each reads quality information
        q_score_list = Vector{Number}()
        #Store dimer information
        dimer_count = 0

        fastq = GZip.open(fastq_file)
        sample_read_count = 1
        line_tracker = 1

        for line in eachline(fastq)
            #Ignore lines without sequence or quality information
            if line_tracker % 4 != 2 && line_tracker % 4 != 0
                line_tracker += 1
            elseif line_tracker % 4 == 2
                #Count canonical 0 bp dimer
                dimer = startswith(line, "TGGAATTCTCGGGTGCC")

                if dimer
                    dimer_count += 1
                end

                sample_read_count += 1
                line_tracker += 1
            elseif line_tracker % 4 == 0
                #Calculate average quality score
                GetReadQScore!(line, q_score_list)
                line_tracker += 1
            end  
        end

        close(fastq)

        #Add read count to dictionary for use later
        read_count_dict[sample_name] = sample_read_count

        #Calculate q score for each sample
        average_q_score = round(mean(q_score_list), digits = 2)
        q_score_dict[sample_name] = average_q_score

        percent_dimer = 100 * dimer_count / sample_read_count
        dimer_count_dict[sample_name] = percent_dimer

        next!(progress_bar_update)
    end

    return read_count_dict, dimer_count_dict, q_score_dict
end

"""
    TrimAdapters(Fastqs::Vector{String}
                    , Sample_Names::Vector{SubString{String}}
                    )

Trim the 3' adapter from each read.

Combo-Seq Read 1 Setup:\n
               5' Adapter      -          Insert        - 3' Adapter      
GATCGTCGGACTGTAGAACTCTGAACNNNN - TGTCAGTTTGTCAAATACCCCA - AAAAAAAAAA 

The bases on the 3' end are also quality trimmed if their quality score is below 20. Reads 
shorter than 16 bases or that weren't trimmed are discarded.

#Example
```julia
julia> TrimAdapters(["sub_sample1.fastq","sub_sample2.fastq","sub_sample3.fastq"]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1.cut.fastq"
 "sample2.cut.fastq"
 "sample3.cut.fastq"
```
"""
function TrimAdapters(Fastqs::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )
    number_of_records = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Trimming adapters...")

    for (fastq_file, sample_name) in zip(Fastqs, Sample_Names)
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

        next!(progress_bar_update)
    end

    trimmed_fastq_files::Vector{String} = CaptureTargetFiles(".cut.fastq")

    return trimmed_fastq_files
end

"""
    AlignWithSalmon(Trimmed_Fastq_Files::Vector{String}
                        , Salmon_Reference::AbstractString
                        , Sample_Names::Vector{SubString{String}}
                        , Need_Reference::AbstractString
                        )
This function aligns trimmed FASTQ files to a reference genome using Salmon.

The library type (libType) is set to stranded forward (SF) since the Combo-Seq protocol is
stranded in that direction.

The function takes the following arguments:

* `Trimmed_Fastq_Files`: A vector of strings containing the paths to the trimmed FASTQ files.
* `Salmon_Reference`: A string containing the path to the Salmon reference genome.
* `Sample_Names`: A vector of strings containing the names of the samples.
* `Need_Reference`: A string indicating whether or not a custom Salmon index was created.

The function returns a vector of strings containing the paths to the Salmon output files.
"""
function AlignWithSalmon(Trimmed_Fastq_Files::Vector{String}
                        , Salmon_Reference::AbstractString
                        , Sample_Names::Vector{SubString{String}}
                        , Need_Reference::AbstractString
                        )
    number_of_records::Int64 = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Aligning reads with Salmon...")

    #Iterate over the fastq files
    for (fastq_file, sample_name) in zip(Trimmed_Fastq_Files, Sample_Names)

        #Check if a custom index was created
        if Need_Reference == "y"
            salmon_reference = string("data/", Salmon_Reference)
        else
            salmon_reference = Salmon_Reference
        end

        #Align the fastq file to the salmon reference
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

        #Rename the output files
        try
            mv("quant.sf", "$sample_name.quant.sf", force=true)
            mv("logs/salmon_quant.log", "$sample_name.salmon_quant.log", force=true)
            mv("lib_format_counts.json", "$sample_name.lib_format_counts.json", force=true)
        catch
            println("Error, no Salmon output found.")
            println("Exiting program...")
            exit()
        end

        #Remove the transcript aligned sam files
        TrashRemoval(CaptureTargetFiles(".mRNA.sam"))

        #Update the progress bar
        next!(progress_bar_update)
    end

    #Get the names of the salmon output files
    salmon_mrna_counts = CaptureTargetFiles("quant.sf")

    return salmon_mrna_counts
end

"""
Use aligned SAM files to create a read length distribution file.
"""
function CalculateReadLengthDistribution(SAM_Files::Vector{String}
                                        , Sample_Names::Vector{SubString{String}}
                                        )
    number_of_records = length(SAM_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating Read Length Distribution..."
                                            )

    for (sam_file, sample_name) in zip(SAM_Files, Sample_Names)
	    #=
        If the SAM file is empty or almost empty, it will cause errors 
        downstream so the script exits.
        =#
        read_length_distibution = 
        try
            read(pipeline(
                `samtools 
                stats 
                -@ 12
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

        next!(progress_bar_update)

        #Write read length data to output file
        output_length_file::IOStream = open("read_lengths_$sample_name.tsv", "w")
        write(output_length_file, string("Length", "\t", "Reads", "\n"))
        write(output_length_file, read_length_distibution)

        close(output_length_file)
    end

    length_files::Vector{String} = CaptureTargetFiles("read_lengths_")

    return length_files
end

"""
Create barplot of each sample's fragment lengths based on the number of reads found at each length
"""
function PlotFragmentLengths(Length_Files::Vector{String}
                            , Sample_Names::Vector{SubString{String}}
                            )
    number_of_records = length(Length_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , 1
                                            , "Creating Read Length Plot..."
                                            )

    for (file, sample_name) in zip(Length_Files, Sample_Names)
        length_file::DataFrame = CSV.read(file, DataFrame)

        #Isolate fragment lengths and the number of reads at each length
        fragment_lengths::Vector{Int64} = length_file[!, :Length]
        reads_per_length::Vector{Int64} = length_file[!, :Reads]
        next!(progress_bar_update)

        plot = Figure(resolution = (1800, 1200))
        next!(progress_bar_update)

        Axis(plot[1, 1]
            , xticks = 15:5:maximum(fragment_lengths)
            , xlabel = "Fragment Lengths"
            , ylabel = "Number of Reads"
            , title = sample_name
            )
        next!(progress_bar_update)

        #Make sample barplot
        barplot!(fragment_lengths
                , reads_per_length
                , fillto = -1
                , color = fragment_lengths
                , strokecolor = :black
                , strokewidth = 1
                )

        save(string(sample_name, "_fragment_lengths.pdf"), plot)
        next!(progress_bar_update)
    end

    pdf_length_files::Vector{String} = CaptureTargetFiles("_fragment_lengths.pdf")

    merge_pdfs([pdf_length_files...], "ComboSeq_fragment_length_plots.pdf")

    TrashRemoval(pdf_length_files)
end

function CalculateDuplicateReads(Fastq::String)
    read_counter = 0
    fastq::IOStream = open(Fastq, "r")
    all_sequences = Vector{String}()

    for line in eachline(fastq)
        #Add every second line from the fastq, .i.e. the read sequence
        push!(all_sequences, readline(fastq))

        #Skip spacer and quality strings
        readline(fastq)   
        readline(fastq)

        #Track the number of reads in the fastq
        read_counter += 1
    end

    close(fastq)

    unique_reads = 100 * length(unique(all_sequences)) / read_counter

    return round(100 - unique_reads, digits = 2)
end

"""
    CalculateSalmonMetrics(Trimmed_Fastq_Files::Vector{String}
                            , Read_Count_Dict::Dict{String, Int64}
                            , Sample_Names::Vector{SubString{String}}
                            , Dimer_Count_Dict::Dict{String, Float64}
                            , Q_Score_Dict::Dict{String, Number}
                            )

Align fastq with specified reference RNA types.

Divide reads aligned with total reads to calculate percent alignment. Store 
alignment information in a dictionary.

#Example
```julia
julia> CalculateSalmonMetrics("sample1.cut.fastq"
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928))
Dict{String, Float64} with 5 entries:
  "miRNA" => 65.0
  "tRNA" => 4.2
  "piRNA" => 1.2
  "snoRNA" => 1.9
  "rRNA" => 3.3
```
"""
function CalculateSalmonMetrics(Trimmed_Fastq_Files::Vector{String}
                                , Read_Count_Dict::Dict{String, Int64}
                                , Sample_Names::Vector{SubString{String}}
                                , Dimer_Count_Dict::Dict{String, Float64}
                                , Q_Score_Dict::Dict{String, Number}
                                )
    number_of_records = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating Sample Metrics..."
                                            )

    sample_tracker = 1

    #Make output metrics file to write stats into
    metrics_file = open("ComboSeq_metrics.csv", "w")

    #Write metrics header information
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
        , "Duplicates", ","
        , "Unaligned Reads"
        , "\n"
        ))

    for (fastq_file, sample_name) in zip(Trimmed_Fastq_Files, Sample_Names)
        #Calculate q_score score
        average_q_score = Q_Score_Dict[sample_name]

        #Gather the number of canonical dimer reads
        percent_dimer = Dimer_Count_Dict[sample_name]

        #Calculate the percent duplicates in each sample
        deduplication_statistics = CalculateDuplicateReads(fastq_file)

        #Gather metrics data from accumlulated analysis files
        cutadapt_information = read(`cat $sample_name.cutadapt_information`, String)
        salmon_library_counts_file = read(`cat $sample_name.lib_format_counts.json`, String)
        salmon_log_file = read(`cat $sample_name.salmon_quant.log`, String)
        miRNA_aligned_reads = readchomp(`samtools view -c -F 0x4 $sample_name.miRNA.sam`)

        #Uses Salmon output JSON file which reports the number of fragments that had 
        #at least one mapping compatible with a stranded library
        total_directional_RNAs_string = match(r"\"num_assigned_fragments\":.\K\d+", salmon_library_counts_file).match
        consistent_directional_RNAs_string = match(r"\"num_frags_with_concordant.+\":.\K\d+", salmon_library_counts_file).match
        
        #Uses Salmon output log file to gather the mRNA mapping rate
        mRNA_mapping_string = match(r"Mapping rate =.\K\d+\.\d{2}", salmon_log_file).match
        
        #Search cutadapt trimming information so that adapter and short fragment data can be easily captured by Regex
        short_fragments_string = match(r"Reads that were too short:.+\(\K.+(?=\%)", cutadapt_information).match
        
        #Convert gathered strings into floats for calculations and adding to the output file
        #They need to be floats in order to plot them later
        miRNA_aligned_reads = parse(Float64, miRNA_aligned_reads)
        short_fragments = abs(parse(Float64, short_fragments_string) - percent_dimer)
        mRNA_mapping = parse(Float64, mRNA_mapping_string)
        total_directional_RNAs = parse(Float64, total_directional_RNAs_string)
        consistent_directional_RNAs = parse(Float64, consistent_directional_RNAs_string)

        #Calculate strand directionality
        percent_directionality::Float64 = round(100 * consistent_directional_RNAs / total_directional_RNAs, digits = 2)

        #Calculate miRNA alignment rate
        miRNA_mapping::Float64 = round(100 * miRNA_aligned_reads / Read_Count_Dict[sample_name], digits = 2)

        aligned_reads = round(miRNA_mapping + mRNA_mapping, digits = 2)
        unaligned_reads = round(100 - aligned_reads - percent_dimer - short_fragments, digits = 2)

        #Write metrics data to output file
        write(metrics_file, string(
            sample_name, ","
            , Read_Count_Dict[sample_name], ","
            , aligned_reads, ","
            , average_q_score, ","
            , round(percent_dimer, digits = 2), ","
            , round(short_fragments, digits = 2), ","
            , miRNA_mapping, ","
            , mRNA_mapping, ","
            , percent_directionality, ","
            , deduplication_statistics, ","
            , unaligned_reads
            , "\n"
            ))

        next!(progress_bar_update)
    end

    close(metrics_file)

    #Returns the metrics file name
    return "ComboSeq_metrics.csv"
end

"""
Create violin plot of each sample's RNA aligment metrics
"""
function MakeMetricsViolinPlot(Metrics_File::String)

    #Import sample metrics file
    metrics_file::DataFrame = CSV.read(Metrics_File, DataFrame)

    #Create new dataframe without sample names and add sample and column names to arrays
    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    number_of_records = length(column_names)
    progress_bar_update = ProgressBarUpdate(number_of_records, .25
                                            , "Creating Violin Plots..."
                                            )
    
    #Gather each metric's column data for plotting
    metrics_columns(a) = [number for number in metrics_file[!, a]]

    #Make empty figure and set x and y axes
    #x-axis is vector of integers from 1 to the number of metrics
    #y-axis is vector of vectors containing the data from each metric
    fig = Figure(resolution = (1800, 1200))
    x, y = 1:num_of_cols, [metrics_columns(column) for column in 1:num_of_cols]
    
    #Create violin plot with data from all samples analyzed
    for column in 1:num_of_cols
        ax = Axis(fig[1, column]
                , yticks = 0:5:100
                , ylabel = "Percent"
                , title = column_names[column]
                )
        CairoMakie.ylims!(ax, 0, 100)

        #Make violin plot with combined sample data
        CairoMakie.violin!(fig[1, column]
                            , repeat([x[column]]
                            , first(size(metrics_file))
                            )
                            , y[column]
                            , show_median=true
                            )
        next!(progress_bar_update)

        #Save file once all columns have been added
        if column == num_of_cols
            save("Violin_plot_metrics.png", fig)
        end

        next!(progress_bar_update)
    end
end

"""
Create barplot of each sample's RNA aligment metrics.
"""
function PlotMetrics(Metrics_File::String
                    , Sample_Names::Vector{SubString{String}}
                    )

    metrics_file::DataFrame = CSV.read(Metrics_File, DataFrame)

    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    colors = [:snow3, :honeydew2, :lightpink, :bisque2, :deepskyblue, :aquamarine3
            , :peachpuff, :paleturquoise1, :orange1, :lightsalmon2, :lemonchiffon
            ]
    
    sample_rows(a) = [values for values in metrics_file[a, :]]

    number_of_records = length(Sample_Names)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Creating Metrics Plots..."
                                            )
    
    for sample in eachindex(Sample_Names)
        sample_name = Sample_Names[sample]

        #Set axis values; y-axis is a vector containing each sample's calculated metrics
        x, y = 1:num_of_cols, sample_rows(sample)
        fig = Figure(resolution = (1800, 1200))

        ax = Axis(fig[1, 1]
                , xticks = (1:num_of_cols
                , column_names)
                , yticks = 0:10:100
                , ylabel = "Percent"
                , title = string(sample_name, " Metrics - ", "ComboSeq")
                )

        #Set y-axis limits
        CairoMakie.ylims!(ax, 0, 105)

        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , strokecolor = :black
                , color = colors[1:num_of_cols]
                , bar_labels = :y
                )

        next!(progress_bar_update)

        save(string(sample_name, "_metrics.pdf"), fig)
    end

    pdf_metrics_files::Vector{String} = CaptureTargetFiles("_metrics.pdf")

    merge_pdfs([pdf_metrics_files...], "ComboSeq_metrics_plots.pdf")
    
    TrashRemoval(pdf_metrics_files)
end

"""
Create a file with a list of unique species names given a miRBase mature or hairpin fasta.
"""
function GetGenusNames(miRBase_Fasta::String)
    miRNA_reference = GZip.open(miRBase_Fasta)
    genus_dictionary = Dict{String, String}()

    for line in eachline(miRNA_reference)
        genus_abbreviation = match(r">\K\w+", line)
        genus_name = match(r">.+MI\w+\s\K\w+", line)
        if !isnothing(genus_name)
            genus_dictionary[genus_name.match] = genus_abbreviation.match
        end
    end

    close(miRNA_reference)

    open("mirBase_genus_list.txt", "w") do output_genus_file
		for (genus, abbreviation) in genus_dictionary
			write(output_genus_file, string(genus, ",", abbreviation, "\n"))
		end
	end

end

"""
Function takes in mirBase mature miRNA fasta file with data from all available
organisms and pulls out only the miRNA data pertaining to the target organism.
"""
function CreateTargetOrganismFastaFile(Organism_Name::String)

    #Skip making species miRNA reference if it exists
    bowtie2_reference_files = CaptureTargetFiles(".bt2")
    if isempty(bowtie2_reference_files)
        #Download miRBase mature miRNA sequence file
        wait(run(pipeline(`wget 
                    https://www.mirbase.org/ftp/CURRENT/mature.fa.gz 
                    --no-check-certificate`
                    , devnull)
                    , wait = false
                )
            )

        #Create file with all genuses with miRBase annotations
        GetGenusNames("mature.fa.gz")

        #Sets the target or organism's genus and species for labeling the output file
        organism_genus_name::String = first(split(Organism_Name, " "))
        organism_species_name::String = last(split(Organism_Name, " "))
        target_genus_abbreviation = ""
        output_fasta_file_name = organism_genus_name * "_" * organism_species_name * ".fa"
        bowtie2_reference_name = organism_genus_name * "_" * organism_species_name

        #Opens the IO streams and sets the output file names for the fasta and bowtie2 index
        genus_file = open("mirBase_genus_list.txt", "r")
        input_fasta_file = GZip.open("mature.fa.gz")
        output_fasta_file = open(organism_genus_name * "_" * organism_species_name * ".fa", "w")

        #Loops through file containing the genus of all the organisms with miRNA data
        #If the target organism is in that list, the genus abbreviation is stored for use below
        for line in eachline(genus_file)
            if occursin(lowercase(organism_genus_name), lowercase(line))
                target_genus_abbreviation = last(split(line, ","))
                break
            end
        end

        #=
        Loops through mirBase mature miRNA fasta looking for the header and sequence information
        for the target genus. If there's a match, that information is added to a new fasta file.
        Must convert RNA sequences in miRBase fasta to DNA.
        =#
        for line in eachline(input_fasta_file)
            input_genus_abbreviation = match(r">\K\w+", line)
            if !isnothing(input_genus_abbreviation)
                if target_genus_abbreviation == input_genus_abbreviation.match
                    write(output_fasta_file, line, "\n")
                    write(output_fasta_file, replace(readline(input_fasta_file), "U" => "T"))
                end
            end
        end

        close(input_fasta_file)
        close(output_fasta_file)
        close(genus_file)

        #=
        Uses the newly created single organism fasta to build a bowtie2 index; will be used
        for alignment with bowtie2 to determine miRNA count information.
        =#
        wait(run(pipeline(
            `bowtie2-build 
            $output_fasta_file_name 
            $bowtie2_reference_name`
            , devnull)
            , wait = false
            )
        )

        TrashRemoval("mirBase_genus_list.txt")
        TrashRemoval("mature.fa.gz")
    end
end

"""
    GenerateMiRNACounts(SAM_File::String)

This function creates a DataFrame containing all aligned miRNA and the number of times they appear.

The function takes the following arguments:

* `SAM_File`: A string containing the path to the SAM file.

The function returns a DataFrame with the following columns:

* `name`: The name of the miRNA.
* `count`: The number of times the miRNA appears in the SAM file.
"""
function GenerateMiRNACounts(SAM_File::String)
    miRNA_names = Vector{String}()
    #Open the SAM file and read each line
    open(SAM_File, "r") do sam_file
        for line in eachline(sam_file)
            miRNA_name = match(r"0\s+\K[a-z][a-z-A-Z0-9]+", line)
            #Check if the line contains a miRNA name
            if !isnothing(miRNA_name)
                #Add the miRNA name to the list of miRNA names
                push!(miRNA_names, miRNA_name.match)
            end
        end
    end

    #Count the number of times each miRNA name appears
    miRNA_counts = countmap(miRNA_names)

    #Create a DataFrame with the miRNA names and counts
    miRNA_counts_dataframe = DataFrame([collect(keys(miRNA_counts))
                                        , collect(values(miRNA_counts))]
                                        , [:name, :count]
                                        )

    #Sort the DataFrame by count in descending order
    sort!(miRNA_counts_dataframe, :2, rev = true)

    return miRNA_counts_dataframe
end

"""
    miRNADiscoveryCalculation(Trimmed_Fastq_Files::Vector{String}
                                , Sample_Names::Vector{SubString{String}}
                                , Bowtie2_Reference_Name::AbstractString
                                , Organism_Name::AbstractString
                                , Need_Reference::AbstractString
                                )

Align trimmed fastq files to the single organism miRNA bowtie2 reference. 

The bowtie2 output SAM is input into function to generate miRNA read counts.

#Example
```julia
julia> miRNADiscoveryCalculation(["sample1.cut.fastq","sample2.cut.fastq","sample3.cut.fastq"]
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
function miRNADiscoveryCalculation(Trimmed_Fastq_Files::Vector{String}
                                    , Sample_Names::Vector{SubString{String}}
                                    , Bowtie2_Reference_Name::AbstractString
                                    , Organism_Name::AbstractString
                                    , Need_Reference::AbstractString
                                    )
    number_of_records = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating the number of miRNA present..."
                                            )

    vector_of_miRNA_counts_dfs = Vector{DataFrame}()

    for (fastq_file, sample_name) in zip(Trimmed_Fastq_Files, Sample_Names)

        if Need_Reference == "y"
            bowtie2_reference = string("data/", Bowtie2_Reference_Name)
        else
            bowtie2_reference = string("data/", join(split(Organism_Name), "_"))
        end

        #Align to miRBase reference
        wait(run(pipeline(
                `bowtie2
                --norc
                --threads 12
                -x $bowtie2_reference
                -U $fastq_file
                -S $sample_name.miRNA.sam`
                , devnull)
                , wait = false
                )
            )

        counts_df = GenerateMiRNACounts(string(sample_name, ".miRNA.sam"))
        push!(vector_of_miRNA_counts_dfs, counts_df)        
        next!(progress_bar_update)
    end

    sam_files::Vector{String} = CaptureTargetFiles(".miRNA.sam")

    return vector_of_miRNA_counts_dfs, sam_files
end

"""
    PlotMiRNACounts(miRNA_Counts_Dfs::Vector{DataFrame}
                        , Sample_Names::Vector{SubString{String}}
                        )

Plot counts of the number of unique miRNA each sample aligned to.

#Example
```julia
julia> PlotMiRNACounts([sample1_counts_df,sample2_counts_df,sample3_counts_df]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1_miRNA_counts.csv"
 "sample2_miRNA_counts.csv"
 "sample3_miRNA_counts.csv"
```
"""
function PlotMiRNACounts(miRNA_Counts_Dfs::Vector{DataFrame}
                        , Sample_Names::Vector{SubString{String}}
                        )
    number_of_records = length(miRNA_Counts_Dfs)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Counting the miRNA in each sample..."
                                            )

    for (miRNA_df, sample_name) in zip(miRNA_Counts_Dfs, Sample_Names)
        #Calculate the total number of reads mapped to miRNA
        total_mapped_reads = sum(miRNA_df[!, :count])

        #Determine reads per million (RPM) for each sample
        RPM::DataFrame = combine(miRNA_df
        , :count => ByRow(miRNA_reads -> round(miRNA_reads / total_mapped_reads * 10^6, digits = 2)
        ))

        #Remove "(unk)" from all miRNA names
        clean_miRNA_names::DataFrame = combine(miRNA_df
                                                , :name => ByRow(name -> replace(name, "(unk)" => ""))
                                                )

        #Reconstitute dataframe with new miRNA names, miRNA counts, and RPM column
        miRNA_df::DataFrame = hcat(clean_miRNA_names, select(miRNA_df, :count), RPM)

        #Rename column with miRNA names and RPM data
        rename!(miRNA_df, :count_function => :RPM, :name_function => :name)

        #Check column containing miRNA counts and adds miRNA above a set threshold
        threshold_of_one = length(filter(>=(1), miRNA_df[!, :count]))
        threshold_of_three = length(filter(>=(3), miRNA_df[!, :count]))
        threshold_of_five = length(filter(>=(5), miRNA_df[!, :count]))
        threshold_of_ten = length(filter(>=(10), miRNA_df[!, :count]))

        #Barplot colors
        colors = [:grey88, :skyblue2, :peachpuff, :lightsalmon]

        #Set x and y-axis values; y-axis is a vector containing each sample's counted miRNAs
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

        #Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
                )

        save(string(sample_name, "_miRNA_counts.pdf"), fig)

        #Write output file containing miRNA names, read count, and RPM data
        CSV.write(string(sample_name, "_miRNA_counts.csv"), miRNA_df)
        
        next!(progress_bar_update)
    end

    pdf_miRNA_files::Vector{String} = CaptureTargetFiles("_miRNA_counts.pdf")

    #Merge output PDFs into one file
    merge_pdfs([pdf_miRNA_files...], "Comboseq_miRNA_plots.pdf")
    
    #Remove individual sample plots
    TrashRemoval(pdf_miRNA_files)

    miRNA_counts_files::Vector{String} = CaptureTargetFiles("_miRNA_counts.csv")

    return miRNA_counts_files
end

function PlotMRNACounts(mRNA_Files::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )
    number_of_records = length(mRNA_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Counting the mRNA in each sample..."
                                            )

    for (file, sample_name) in zip(mRNA_Files, Sample_Names)
        #Import Salmon output file containing mRNA information
        mRNA_file::DataFrame = CSV.read(file, DataFrame)

        #Checks column containing mRNA counts and adds mRNA above a set threshold
        threshold_of_one = sum(mRNA_file[!,:TPM] .>= 1)
        threshold_of_ten = sum(mRNA_file[!,:TPM] .>= 10)
        threshold_of_one_hundred = sum(mRNA_file[!,:TPM] .>= 100)
        threshold_of_five_hundred = sum(mRNA_file[!,:TPM] .>= 500)
        threshold_of_one_thousand = sum(mRNA_file[!,:TPM] .>= 1000)

        #Captures the name and TPM count of each mRNA with a TPM of at least 1
        mRNA_at_threshold_of_one = mRNA_file[mRNA_file[!,:TPM] .>= 1, [:Name, :TPM]]

        #Sort output dataframe with mRNA TPM counts by the most highly expressed mRNAs
        sorted_mRNA_at_threshold_of_one = sort(mRNA_at_threshold_of_one[!, [:Name, :TPM]]
                                                , :TPM
                                                , rev = true
                                                )

        sample_name = Sample_Names[sample_tracker]

        #Barplot colors
        colors = [:burlywood, :lightsteelblue, :mistyrose2, :pink4, :grey80]

        #Set x and y-axis values; y-axis is a vector containing each sample's counted mRNAs 
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

        #Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
                )

        save(sample_name * "_mRNA_counts.pdf", fig)

        #Write output file containing mRNA names and TPM data
        CSV.write(sample_name * "_mRNA_counts.csv", sorted_mRNA_at_threshold_of_one)
        
        next!(progress_bar_update)
    end

    pdf_mRNA_files::Vector{String} = CaptureTargetFiles("_mRNA_counts.pdf")

    #Merge output PDFs into one file
    merge_pdfs([pdf_mRNA_files...], "Comboseq_mRNA_plots.pdf")
    
    #Remove individual sample plots
    TrashRemoval(pdf_mRNA_files)

    mRNA_counts_files::Vector{String} = CaptureTargetFiles("_mRNA_counts.csv")

    return mRNA_counts_files
end

"""
Find all miRNAs the samples have in common. Calculate the reads per million reads (RPM).
"""
function FindCommonRNAs(miRNA_Counts_Files::Vector{String}
                        , Read_Count_Dict::Dict{String, Int64}
                        , Sample_Names::Vector{SubString{String}}
                        )
    #Vector of dictionaries containing the miRNA counts from each sample
    miRNA_info = Vector{Dict{String, Int64}}()
    RPM_info = Vector{Dict{String, Number}}()
    full_miRNA_names_list = Vector{String}()

    for (file, sample_name) in zip(miRNA_Counts_Files, Sample_Names)
        miRNA_file::IOStream = open(file, "r")

        #Dictionaries to hold each miRNA and its associated read count
        sample_miRNA_counts_dictionary = Dict{String, Int64}()
        sample_miRNA_RPM_dictionary = Dict{String, Number}()

        #Skip header line
        readline(miRNA_file)

        for line in eachline(miRNA_file)
            split_line = split(line, ",")
            miRNA_count = parse(Int64, split_line[2])
            miRNA_name = first(split_line)
            RPM = round(miRNA_count / (Read_Count_Dict[sample_name] / 10^6), digits = 4)

            sample_miRNA_counts_dictionary[miRNA_name] = miRNA_count
            sample_miRNA_RPM_dictionary[miRNA_name] = RPM
        end

        push!(miRNA_info, sample_miRNA_counts_dictionary)
        push!(RPM_info, sample_miRNA_RPM_dictionary)
        
        append!(full_miRNA_names_list, keys(sample_miRNA_counts_dictionary))

        close(miRNA_file)
    end

    #Count number each miRNA's occurances
    miRNA_names_dict = countmap(full_miRNA_names_list)

    #Find miRNA present in all samples, .i.e. have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(Sample_Names), miRNA_names_dict) |> keys

    return miRNAs_in_common, miRNA_info, RPM_info
end

"""
Find all RNAs the samples have in common.
"""
function FindCommonRNAs(mRNA_Counts_Files::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )
    #Vector of dictionaries containing the mRNA counts from each sample
    mRNA_info = Vector{Dict{String, Float64}}()
    TPM_info = Vector{Dict{String, Number}}()
    full_mRNA_names_list = Vector{String}()

    for file in mRNA_Counts_Files
        mRNA_file::IOStream = open(file, "r")

        #Dictionaries to hold each mRNA and its associated read count
        sample_mRNA_counts_dictionary = Dict{String, Float64}()
        sample_mRNA_TPM_dictionary = Dict{String, Number}()

        #Skip header line
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

    #Count number each mRNA's occurances
    miRNA_names_dict = countmap(full_mRNA_names_list)

    #Find miRNA present in all samples, .i.e. have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(Sample_Names), miRNA_names_dict) |> keys

    return miRNAs_in_common, mRNA_info, TPM_info
end

function WriteCommonRNAFile(miRNA_Names::Base.KeySet{String, Dict{String, Int64}}
                            , miRNA_Info::Vector{Dict{String, Int64}}
                            , RPM_Info::Vector{Dict{String, Number}}
                            , Sample_Names::Vector{SubString{String}}
                            ; RNA_Type="miRNA"
                            )

    for (index, sample) in enumerate(Sample_Names)
        common_miRNA_file::IOStream = open(string(sample, "_common_miRNAs.tsv"), "w")
        common_miRNA_file_RPM::IOStream = open(string(sample, "_common_miRNAs_RPM.tsv"), "w")

        if index == 1
            write(common_miRNA_file, string(RNA_Type, "\t", sample, "\n"))
            write(common_miRNA_file_RPM, string(RNA_Type, "\t", sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(miRNA, "\t", miRNA_Info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(miRNA, "\t", RPM_Info[index][miRNA], "\n"))
            end
        else
            write(common_miRNA_file, string(sample, "\n"))
            write(common_miRNA_file_RPM, string(sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(miRNA_Info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(RPM_Info[index][miRNA], "\n"))
            end
        end

        close(common_miRNA_file)
        close(common_miRNA_file_RPM)
    end

    common_files = CaptureTargetFiles("_common_miRNAs.tsv")
    common_RPM_files = CaptureTargetFiles("_common_miRNAs_RPM.tsv")

    #Combine all the separate common miRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_miRNAs.tsv"))
    run(pipeline(`paste $common_RPM_files`, stdout = "Common_RPM_miRNAs.tsv"))

    #Run principal component analysis and UMAP on common miRNA
    PlotClustering("Common_miRNAs.tsv", RNA_Type)

    TrashRemoval(common_files)
    TrashRemoval(common_RPM_files)
end

function WriteCommonRNAFile(mRNA_Names::Base.KeySet{String, Dict{String, Int64}}
                            , mRNA_Info::Vector{Dict{String, Float64}}
                            , TPM_Info::Vector{Dict{String, Number}}
                            , Sample_Names::Vector{SubString{String}}
                            ; RNA_Type="mRNA"
                            )
    
    for (index, sample) in enumerate(Sample_Names)
        common_mRNA_file::IOStream = open(string(sample, "_common_mRNAs.tsv"), "w")
        common_mRNA_file_TPM::IOStream = open(string(sample, "_common_mRNAs_TPM.tsv"), "w")

        if index == 1
            write(common_mRNA_file, string(RNA_Type, "\t", sample, "\n"))
            write(common_mRNA_file_TPM, string(RNA_Type, "\t", sample, "\n"))
            for mRNA in mRNA_Names
                write(common_mRNA_file, string(mRNA, "\t", mRNA_Info[index][mRNA], "\n"))
                write(common_mRNA_file_TPM, string(mRNA, "\t", TPM_Info[index][mRNA], "\n"))
            end
        else
            write(common_mRNA_file, string(sample, "\n"))
            write(common_mRNA_file_TPM, string(sample, "\n"))
            for mRNA in mRNA_Names
                write(common_mRNA_file, string(mRNA_Info[index][mRNA], "\n"))
                write(common_mRNA_file_TPM, string(TPM_Info[index][mRNA], "\n"))
            end
        end

        close(common_mRNA_file)
        close(common_mRNA_file_TPM)
    end

    common_files = CaptureTargetFiles("_common_mRNAs.tsv")
    common_TPM_files = CaptureTargetFiles("_common_mRNAs_TPM.tsv")

    #Combine all the separate common mRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_mRNAs.tsv"))
    run(pipeline(`paste $common_TPM_files`, stdout = "Common_TPM_mRNAs.tsv"))

    #Run principal component analysis and UMAP on common mRNA
    PlotClustering("Common_mRNAs.tsv", RNA_Type)

    TrashRemoval(common_files)
    TrashRemoval(common_TPM_files)
end

"""
Plot principal component analysis (PCA) and RNA counts UMAP of all libraries if possible.
"""
function PlotClustering(Common_RNA_File::String, RNA_Type::String)
    common_RNA_counts = DataFrame(CSV.File(Common_RNA_File))
    sample_names = DataFrames.names(common_RNA_counts, Not([Symbol("$RNA_Type")]))
    RNA_names = common_RNA_counts[!, Symbol("$RNA_Type")]

    #There must be at least two samples in order to perform the principal component analysis
    #Must also have at least one RNA in common
    if last(size(common_RNA_counts)) > 2 && first(size(common_RNA_counts)) > 1
        common_RNA_counts_matrix = Matrix{Float64}(select(common_RNA_counts, Not([Symbol("$RNA_Type")])))

        #Make PCA covariance matrix for PCA plot
        PCA_matrix = fit(PCA, common_RNA_counts_matrix; maxoutdim = 20)
        transformed_counts = predict(PCA_matrix, common_RNA_counts_matrix)

        #Make PCA covariance matrix for UMAP plot
        UMAP_PCA_matrix = fit(PCA, common_RNA_counts_matrix'; maxoutdim = 20)
        transposed_transformed_counts = predict(UMAP_PCA_matrix, common_RNA_counts_matrix')

        #Plot PCA only if there are at least two principal components
        if first(size(transformed_counts)) >= 2
            sample_name_df = DataFrame("samples" => sample_names)
            RNA_names_df = DataFrame("$RNA_Type" => RNA_names)

            #Write PCA information to file for cluster tracking
            PCA_values = hcat(sample_name_df, DataFrame(transformed_counts[1:2, :]', ["PC1", "PC2"]))
            CSV.write("PCA_information.csv", PCA_values)
            
            #Plot PCA
            Plots.scatter(size = (1200, 800), dpi = 300, titlefont = (16, "Computer Modern")
                            , xlabel = "PC1", ylabel = "PC2", title = "Common $RNA_Type: PCA"
                            , transformed_counts[1, :], transformed_counts[2, :]
                            , left_margin = 23mm, right_margin = 8mm, bottom_margin = 8mm
                            , leg = false
                            )
           savefig(string("common_", RNA_Type, "_PCA.png"))
        end
        
        #Plot UMAP only if there are at least three principal components
        if first(size(transposed_transformed_counts)) >= 3
            sample_name_df = DataFrame("samples" => sample_names)
            RNA_names_df = DataFrame("$RNA_Type" => RNA_names)

            #Create low dimensional embedding for UMAP
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

            #Write UMAP information to file for cluster tracking
            UMAP_values = hcat(RNA_names_df, DataFrame(embedding', ["UMAP1", "UMAP2"]))
            CSV.write("UMAP_information.csv", UMAP_values)

            #=
            Plot UMAP of common RNAs. Color UMAP clusters based on K-medoids clustering; 
            choose number of clusters based on the number of principal components.
            =#
            distance_matrix = pairwise(Euclidean(), embedding, embedding)
            kmedoids_cluster_colors = kmedoids(distance_matrix, first(size(transposed_transformed_counts)))
            Plots.scatter(size = (1200, 800), embedding[1, :], embedding[2, :]
                            , title="Common $RNA_Type: UMAP", left_margin = 13mm
                            , bottom_margin = 10mm, dpi = 300
                            , marker_z = kmedoids_cluster_colors.assignments
                            , color = :lighttest, xlabel = "UMAP1", ylabel = "UMAP2"
                            , titlefont = (16, "Computer Modern"), leg = false
                            , markersize = 9, markerstrokewidth = 0.1
                            )
            savefig(string("common_", RNA_Type, "_UMAP.png"))
        end
    end
end

"""
Spot remove unecessary intermediate files.
"""
function TrashRemoval(Files_to_Delete::Vector{String})
    for file in Files_to_Delete
        rm(file)
    end
end

function TrashRemoval(File_to_Delete::String)
    rm(File_to_Delete)
end

"""
Remove all intermediate files.
"""
function GarbageCollection()
    Files_to_Delete = Set(vcat(
    CaptureTargetFiles("deduplication_statistics")
    ,CaptureTargetFiles(".bt2")
    ,CaptureTargetFiles(".miRNA.")
    ,CaptureTargetFiles("v22_cluster")
    ,CaptureTargetFiles(".sam")
    ,CaptureTargetFiles(".bam")
    ,CaptureTargetFiles("quant")
    ,CaptureTargetFiles("cut")
    ))

    for file in Files_to_Delete
        rm(file)
    end

    json_files = CaptureTargetFiles(".json")
    TrashRemoval(json_files)
end

function julia_main()::Cint

    #Say hello Issac!
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

    print("Would you like to download reference files for the Salmon alignment? (y/n): ")
    need_reference = readline()
    Fastqs = CaptureTargetFiles("_R1_001.fastq.gz")
    Sample_Names = map(sample -> first(split(sample, "_")), Fastqs)
    Salmon_Reference, Organism_Name = CreateReferenceFile(need_reference)
    Read_Count_Dict, Dimer_Count_Dict, Q_Score_Dict = ParseFastqFiles(Fastqs, Sample_Names)
    Trimmed_Fastqs = TrimAdapters(Fastqs, Sample_Names)
    miRNA_Counts_Dfs, SAM_Files = miRNADiscoveryCalculation(Trimmed_Fastqs
                                                            , Sample_Names
                                                            , Salmon_Reference
                                                            , Organism_Name
                                                            , need_reference
                                                            )
    miRNA_Counts_Files = PlotMiRNACounts(miRNA_Counts_Dfs, Sample_Names)
    Salmon_mRNA_Counts = AlignWithSalmon(Trimmed_Fastqs
                                        , Salmon_Reference
                                        , Sample_Names
                                        , need_reference
                                        )
    mRNA_Counts_Files = PlotMRNACounts(Salmon_mRNA_Counts, Sample_Names)
    Length_Files = CalculateReadLengthDistribution(SAM_Files, Sample_Names)
    PlotFragmentLengths(Length_Files, Sample_Names)
    Metrics_File = CalculateSalmonMetrics(Trimmed_Fastqs
                                        , Read_Count_Dict
                                        , Sample_Names
                                        , Dimer_Count_Dict
                                        , Q_Score_Dict
                                        )
    PlotMetrics(Metrics_File, Sample_Names)
    MakeMetricsViolinPlot(Metrics_File)
    Full_miRNA_Names_List, miRNA_Info, RPM_Info = FindCommonRNAs(miRNA_Counts_Files
                                                                , Read_Count_Dict
                                                                , Sample_Names
                                                                )
    WriteCommonRNAFile(Full_miRNA_Names_List, miRNA_Info, RPM_Info, Sample_Names, RNA_Type="miRNA")
    Full_mRNA_Names_List, mRNA_Info, TPM_Info = FindCommonRNAs(mRNA_Counts_Files
                                                                , Sample_Names
                                                                )
    WriteCommonRNAFile(Full_mRNA_Names_List, mRNA_Info, TPM_Info, Sample_Names, RNA_Type="mRNA")
    GarbageCollection()
    println(" ")
    println("Analysis Finished")
    println(" ")
    
    return 0
end

end