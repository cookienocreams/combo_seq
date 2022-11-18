#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
# Combo-Seq pipeline
#
# This script can be used to analyze Combo-Seq libraries using the pseudo-aligner Salmon
# Script inputs: url for desired reference transcriptome and reference genome fasta files
# Script inputs: Name of the species to be analyzed
#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
print("Would you like to download reference files for the Salmon alignment? (y/n): ")
need_reference :: String = readline()

using Downloads
using ProgressMeter
using CairoMakie
using CSV
using DataFrames

#=
Function asks the user if they would like to download a non-human reference
for use by Salmon during the alignment step. If yes, the files are downloaded
and the genome names are extracted to name the Salmon directory.
=#
function CreateReferenceFile(Need_Reference :: String)

    if Need_Reference === "y"
        #Need to gather organism name to see if there are miRNA annotations for use in alignment later 
        print("What is the scientific name of the organism you would like to use?: ")
        organism_name :: String = readline()

        print("Input website url containing the desired transcriptome reference fasta: ")
        reference_transcriptome_fasta :: String = readline()

        #Checks to see if the url contains the secure hypertext protocol, i.e. 'https:'
        #If yes, that is replaced by the non-secure protocol, i.e. http:, to avoid certificate errors
        transcriptome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_transcriptome_fasta)
        if transcriptome_fasta_hypertext_protocol.match == "https:"
            reference_transcriptome_fasta_address = match(r"(?!.+\/\/).+", reference_transcriptome_fasta)
            reference_transcriptome_fasta = "http:" * reference_transcriptome_fasta_address.match
        end

        #Captures file names to be used to refer to the downloaded file
        reference_transcriptome_fasta_name :: String = last(splitpath(reference_transcriptome_fasta))
        unzipped_reference_transcriptome_fasta_name :: String = replace(reference_transcriptome_fasta_name, ".gz" => "")

        #Downloads desired transcriptome fasta file and determines the filename
        Downloads.download(reference_transcriptome_fasta, "$reference_transcriptome_fasta_name")

        print("Input website url containing the desired genomic reference fasta: ")
        reference_fasta :: String = readline()

        genome_fasta_hypertext_protocol = match(r"^.+(?=\/\/)", reference_fasta)
        if genome_fasta_hypertext_protocol.match == "https:"
            reference_genome_fasta_address = match(r"(?!.+\/\/).+", reference_fasta)
            reference_fasta = "http:" * reference_genome_fasta_address.match
        end

        reference_fasta_name :: String = last(splitpath(reference_fasta))
        unzipped_reference_fasta_name :: String = replace(reference_fasta_name, ".gz" => "")
        salmon_reference_name = first(split(reference_fasta_name, "."))

        #Downloads desired genome fasta file and determines the filename
        Downloads.download(reference_fasta, "$reference_fasta_name")

        #Unzips downloaded fasta before running it through Salmon index generation
        run(`gunzip $reference_transcriptome_fasta_name`)
        run(`gunzip $reference_fasta_name`)
        
        #Creates decoy sequences list needed to mask the genome and improve mRNA quantification
        #See Salmon documentation for more info https://salmon.readthedocs.io/en/latest/salmon.html
        genomic_target_names :: String = 
        read(pipeline(
        `cat $unzipped_reference_fasta_name`
        ,`grep "^>"`
        , `cut -d " " -f 1`)
        , String
        )

        #Writes decoy sequences to a new file for use during Salmon index creation
        decoys_output = open("salmon_decoys.txt", "w")
        write(decoys_output, genomic_target_names)
        close(decoys_output)
        run(`sed -i.bak -e 's/>//g' salmon_decoys.txt`)

        #Adds genome fasta to transcriptome fasta for masking during index creation 
        concatenated_reference = read(pipeline(`
        cat $unzipped_reference_transcriptome_fasta_name
        $unzipped_reference_fasta_name
        `))
        concatenated_reference_file = open("gentrome.fa", "w")
        write(concatenated_reference_file, concatenated_reference)
        close(concatenated_reference_file)

        #Creates reference Salmon database for use in further analyses
        for line in eachline(file)
            if occursin("ENST", line) 
                wait(run(pipeline(
                `salmon 
                index 
                -i $salmon_reference_name 
                --transcripts reference_transcriptome_fasta_name 
                -k 21 
                --threads 12 
                --gencode`)
                , wait=false
                ))
                break
            else 
                wait(run(pipeline(
                `salmon 
                index 
                -i $salmon_reference_name 
                --transcripts reference_transcriptome_fasta_name 
                -k 21 
                --threads 12`)
                , wait=false
                ))
                break
            end
        end
        

    else
        #If a custom database isn't used; script defaults to using a human reference
        salmon_reference_name = "gencode.42"
        organism_name = "Homo Sapiens"
    end

    return salmon_reference_name :: String, organism_name, Need_Reference
end

#=
Function list all files in the current directory. It then loops through the directory files
and checks to see if each files contains the desired file's suffix. If it does, that file
is added to an array comprehension for storage.
=#
function CaptureTargetFiles(Files_To_Capture :: String)
    files_in_directory :: Vector{String} = readdir()
    return [file for file in files_in_directory if occursin(Files_To_Capture, file)]
end

#Updates progress bar on the command line when interval is a float
function ProgressBarUpdate(Number_of_Records::Int64, Interval::Float64, Description::String)
    progress_bar_update = Progress(Number_of_Records
    , dt=Interval
    , barglyphs=BarGlyphs("[=> ]")
    , barlen=100
    , desc=Description
    )

    return progress_bar_update
end

#Updates progress bar on the command line when interval is an integer
function ProgressBarUpdate(Number_of_Records::Int64, Interval::Int64, Description::String)
    progress_bar_update = Progress(Number_of_Records
    , dt=Interval
    , barglyphs=BarGlyphs("[=> ]")
    , barlen=100
    , desc=Description
    )

    return progress_bar_update
end

#=
Function downsamples all the raw fastq files to the same number of reads. The value
chosen is the number of reads in the fastq with the fewest reads. Seqtk is then
used to downsample the files.
=#
function DownsampleRawReads(Fastqs :: Vector{String})
    #Sets the number of items to be analyzed for the progress bar
    number_of_records :: Int64 = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Counting reads in each sample...")

    #Counts the number of reads in each input fastq file and sets the target 
    #read count variable with the fewest reads found in the set of fastqs.
    Target_Read_Count :: Int64 = 1000000000
    for fastq_file in Fastqs
        sample_read_count_string :: SubString{String} = readchomp(pipeline(
            `gzip -dc $fastq_file`
            , `awk 'NR % 4 == 2'`
            , `wc -l`
            ))
        #Convert string to integer
        sample_read_count :: Int64 = parse(Int64, sample_read_count_string)
        if sample_read_count < Target_Read_Count
            Target_Read_Count = sample_read_count
        end
        next!(progress_bar_update)
    end

    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Downsampling reads...")

    #Downsamples each fastq file according to the minimum read count calculated above
    for fastq_file in Fastqs
        sample_name :: String = first(split(fastq_file, "_"))
        downsampled_reads :: String = 
        read(
        `seqtk 
        sample -s100 
        $fastq_file 
        $Target_Read_Count`
        , String
        )

        #Writes downsampled fastq to a new file
        downsampled_output= open("sub_$sample_name.fastq", "w")
        write(downsampled_output, downsampled_reads)
        close(downsampled_output)
        next!(progress_bar_update)
    end

    downsampled_fastq_files :: Vector{String} = CaptureTargetFiles("sub_")

    return downsampled_fastq_files
end

#=
Function takes in an array of fastq files. It then loops through those files and trims
the adapters and 4 randomized bases from each read. The bases on the 3' end are also
quality trimmed if their quality score is below 20. Reads shorter than 16 bases or that
weren't trimmed are discarded.
=#
function TrimAdapters(Fastqs :: Vector{String})
    number_of_records :: Int64 = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Trimming adapters...")

    min_length :: Int64 = 16
    for fastq_file in Fastqs
        sample_name :: String = first(splitext(fastq_file))
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
        ))
        next!(progress_bar_update)
    end

    trimmed_fastq_files :: Vector{String} = CaptureTargetFiles(".cut.fastq")

    return trimmed_fastq_files
end

#=
Function takes in an array of adapter trimmed fastq files, the chosen reference name, and
whether or not a custom index was created. It then loops through the trimmed fastq files 
and aligns them to either the chosen reference or the standard human reference using Salmon.
The library type (libType) is set to stranded forward (SF) since the Combo-Seq protocol is
stranded in that direction. The fld flags denote the expected size of the library fragments.
The final Salmon output files are gathered and stored for analysis later.
=#
function AlignWithSalmon(Trimmed_Fastq_Files :: Vector{String}, Salmon_Reference :: String)
    number_of_records :: Int64 = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Aligning reads with Salmon...")

    for fastq_file in Trimmed_Fastq_Files
        sample_name :: String = first(split(fastq_file, "."))
        wait(run(pipeline(
        `salmon 
        quant
        -i $Salmon_Reference
        --libType SF 
        -r $fastq_file
        --validateMappings 
        --fldMean 30
        --fldSD 5
        --writeMappings=$sample_name.sam
        -o ./`
        , devnull)
        , wait=false
        ))

        #Rename files since they are normally given the same name and overwritten each cycle
        mv("quant.sf", "$sample_name.quant.sf")
        mv("logs/salmon_quant.log", "$sample_name.salmon_quant.log")
        mv("lib_format_counts.json", "$sample_name.lib_format_counts.json")
        next!(progress_bar_update)
    end

    #Creates arrays with the Salmon output files for use in generating output data
    salmon_sam_files :: Vector{String} = CaptureTargetFiles(".sam")
    salmon_mRNA_counts :: Vector{String} = CaptureTargetFiles("quant.sf")

    return salmon_sam_files, salmon_mRNA_counts
end

#=
Function takes in Salmon aligned SAM files and converts them to BAM files. The BAM
files are used to create a read length distribution file.
=#
function ConvertSAMToBAM(Salmon_Aligned_SAM_Files :: Vector{String})
    number_of_records :: Int64 = length(Salmon_Aligned_SAM_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Converting SAM files to BAMs...")

    for sam_file in Salmon_Aligned_SAM_Files
        sample_name :: String = first(split(sam_file, "."))
        run(pipeline(
        `samtools
        view
        -@ 12
        -b
        -o $sample_name.bam
        $sam_file`
        , devnull
        ))
        next!(progress_bar_update)
    end

    bam_files :: Vector{String} = CaptureTargetFiles(".bam")

    return bam_files
end

#=
Function takes in Salmon aligned BAM files uses them to create a read length distribution file.
=#
function CalculateReadLengthDistribution(Salmon_BAM_Files :: Vector{String})
    number_of_records :: Int64 = length(Salmon_BAM_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, 1, "Calculating Read Length Distribution...")

    for bam_file in Salmon_BAM_Files
        sample_name :: String = first(split(bam_file, "."))
        read_length_distibution :: String = 
        read(pipeline(
        `samtools 
        stats 
        -@ 12 
        $bam_file`
        , `grep ^RL`
        , `cut -f 2-`)
        , String
        )
        next!(progress_bar_update)

        #Writes read length data to output file
        output_length_file :: IOStream = open("read_lengths_$sample_name.tsv", "w")
        write(output_length_file, "Length" * "\t" * "Reads" * "\n")
        write(output_length_file, read_length_distibution)
        close(output_length_file)
    end

    length_files :: Vector{String} = CaptureTargetFiles("read_lengths_")

    return length_files
end

#Create barplot of each sample's fragment lengths based on the number of reads found at each length
function PlotFragmentLengths(Length_Files :: Vector{String})
    number_of_records :: Int64 = length(Length_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, 1.5, "Creating Read Length Plot...")

    for file in Length_Files
        sample_name = match(r"read_lengths_sub_\K\w+", file).match
        length_file :: DataFrame = CSV.read(file, DataFrame)

        #Isolate fragment lengths and the number of reads at each length
        fragment_lengths :: Vector{Int64} = length_file[!,:Length]
        reads_per_length :: Vector{Int64} = length_file[!,:Reads]
        next!(progress_bar_update)

        plot = Figure(resolution = (1800, 1200))
        next!(progress_bar_update)
        Axis(plot[1, 1]
        , xticks = 20:5:maximum(fragment_lengths)
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

        save(sample_name * "_fragment_lengths.pdf", plot)
        next!(progress_bar_update)
    end

    pdf_length_files :: Vector{String} = CaptureTargetFiles("_fragment_lengths.pdf")

    #Merge output PDFs into one file
    run(`gs 
    -q 
    -dBATCH 
    -dNOPAUSE 
    -sDEVICE=pdfwrite 
    -sOutputFile=ComboSeq_fragment_length_plots.pdf 
    $pdf_length_files`
    )
    #Removes individual sample plots
    TrashRemoval(pdf_length_files)

end

function CalculateDuplicateReads(Fastq :: String)
    read_counter :: Int64 = 0
    fastq :: IOStream = open(Fastq, "r")
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

#Determine basic library metrics and write them to a file for plotting
function CalculateSalmonMetrics(Trimmed_Fastq_Files :: Vector{String})
    number_of_records :: Int64 = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Calculating Sample Metrics...")

    column_names = "Sample\tReads with adapters\tReads too short\tmiRNA\tmRNA\tDirectionality\tDuplicates"

    #Make output metrics file to write stats into
    metrics_file = open("ComboSeq_metrics.tsv", "w")
    write(metrics_file, column_names, "\n")

    for fastq_file in Trimmed_Fastq_Files
        sample_name :: SubString{String} = match(r"sub_\K\w+(?=\.)", fastq_file).match

        #Calculate the percent duplicates in each sample
        deduplication_statistics = CalculateDuplicateReads(fastq_file)

        #Gather metrics data from accumlulated analysis files
        cutadapt_information = read(`cat sub_$sample_name.cutadapt_information`, String)
        Salmon_library_counts_file = read(`cat sub_$sample_name.lib_format_counts.json`, String)
        Salmon_log_file = read(`cat sub_$sample_name.salmon_quant.log`, String)
        miRNA_aligned_reads = readchomp(`samtools view -c -F 0x4 $sample_name.miRNA.bam`)
        total_read_count = readchomp(`samtools view -c $sample_name.miRNA.bam`)

        #Convert gathered strings into floats for calculations and adding to the output file
        #They need to be floats in order to plot them later
        miRNA_aligned_reads = parse(Float64, miRNA_aligned_reads)
        total_read_count = parse(Int64, total_read_count)

        #Calculate miRNA alignment rate
        miRNA_mapping :: Float64 = round(100 * miRNA_aligned_reads / total_read_count, digits = 2)

        #Uses Salmon output JSON file which reports the number of fragments that had at least one mapping compatible with a stranded library
        total_directional_RNAs_string = match(r"\"num_assigned_fragments\":.\K\d+", Salmon_library_counts_file).match
        consistent_directional_RNAs_string = match(r"\"num_frags_with_concordant_consistent_mappings\":.\K\d+", Salmon_library_counts_file).match
        total_directional_RNAs = parse(Float64, total_directional_RNAs_string)
        consistent_directional_RNAs = parse(Float64, consistent_directional_RNAs_string)
        percent_directionality :: Float64 = round(100 * consistent_directional_RNAs / total_directional_RNAs, digits = 2)

        #Uses Salmon output log file to gather the mRNA mapping rate
        mRNA_mapping_string = match(r"Mapping rate =.\K\d+\.\d{2}", Salmon_log_file).match
        mRNA_mapping = parse(Float64, mRNA_mapping_string)

        #Stores cutadapt trimming information in a variable so that adapter and short fragment data can be easily captured by Regex
        reads_with_adapters_string = match(r"Reads with adapters:.+\(\K.+(?=\%)", cutadapt_information).match
        reads_too_short_string = match(r"Reads that were too short:.+\(\K.+(?=\%)", cutadapt_information).match
        reads_with_adapters = parse(Float64, reads_with_adapters_string)
        reads_too_short = parse(Float64, reads_too_short_string)

        #Write data to the output file
        output = "$sample_name\
        \t$reads_with_adapters\
        \t$reads_too_short\
        \t$miRNA_mapping\
        \t$mRNA_mapping\
        \t$percent_directionality\
        \t$deduplication_statistics\
        \n"

        write(metrics_file, output)
        next!(progress_bar_update)

    end

    close(metrics_file)

    #Remove deduplicated fastq
    TrashRemoval("dedup.fastq")

    #Returns the metrics file name
    return match(r"file \K\w+\.\w+", metrics_file.name).match
end

#Create barplot of each sample's RNA aligment metrics
function PlotMetrics(Metrics_File :: String)
    number_of_records :: Int64 = length(Metrics_File) / 2
    progress_bar_update = ProgressBarUpdate(number_of_records, .25, "Creating Metrics Plots...")

    #Import sample metrics file
    full_metrics_file :: DataFrame = CSV.read(Metrics_File, DataFrame)

    #Create new dataframe without sample names and add sample and column names to arrays
    metrics_file :: DataFrame = select(full_metrics_file, Not([:Sample]))
    column_names :: Vector{String} = names(metrics_file)
    sample_names :: Vector{String} = full_metrics_file[!,1]

    #Colors for the bars in the plot
    colors = [:snow3, :honeydew2, :paleturquoise1, :aquamarine3, :peachpuff, :lightsalmon2]
    
    #Function to gather each individual sample row for ploting
    sample_rows(a) = [number for number in metrics_file[a,:]]
    
    for sample in eachindex(sample_names)
        sample_name = sample_names[sample]

        #Set x and y-axis values; y-axis is a vector containing each sample's calculated metrics 
        x, y = 1:size(metrics_file)[2], sample_rows(sample)
        fig = Figure(resolution = (1800, 1200))
        ax = Axis(fig[1, 1]
        , xticks = (1:size(metrics_file)[2]
        , column_names)
        , yticks = 0:10:100
        , ylabel = "Percent"
        , title = sample_name * " Metrics"
        )
        #Set y-axis limits
        ylims!(ax, 0, 105)

        #Make sample barplot
        barplot!(ax
        , x
        , y
        , strokewidth = 1
        , color = colors
        , strokecolor = :black
        , bar_labels = :y
        )
        save(sample_name * "_metrics.pdf", fig)
        next!(progress_bar_update)
    end

    pdf_metrics_files :: Vector{String} = CaptureTargetFiles("_metrics.pdf")

    #Merge output PDFs into one file
    run(`gs 
    -q 
    -dBATCH 
    -dNOPAUSE 
    -sDEVICE=pdfwrite 
    -sOutputFile=ComboSeq_metrics_plots.pdf 
    $pdf_metrics_files`
    )
    #Removes individual sample plots
    TrashRemoval(pdf_metrics_files)

end

#Create violin plot of each sample's RNA aligment metrics
function ViolinPlotMetrics(Metrics_File :: String)

    #Import sample metrics file and delete bottom two rows with labels for the target metrics percentages
    full_metrics_file :: DataFrame = CSV.read(Metrics_File, DataFrame)

    #Create new dataframe without sample names and add sample and column names to arrays
    metrics_file :: DataFrame = select(full_metrics_file, Not([:Sample]))
    column_names :: Vector{String} = names(metrics_file)
    sample_names :: Vector{String} = full_metrics_file[!,1]

    number_of_records :: Int64 = length(column_names)
    progress_bar_update = ProgressBarUpdate(number_of_records, .25, "Creating Violin Plots...")
    
    #Function to gather each metric's column for ploting
    metrics_columns(a) = [number for number in metrics_file[!,a]]

    #Make empty figure and set x and y axes
    #x-axis is vector of integers from 1 to the number of metrics
    #y-axis is vector of vectors containing the data from each metric
    fig = Figure(resolution = (1800, 1200))
    x, y = 1:size(metrics_file)[2], [metrics_columns(column) for column in 1:size(metrics_file)[2]]
    
    #Create violin plot with data from all samples analyzed
    for column in 1:size(metrics_file)[2]
        ax = Axis(fig[1, column]
        , yticks = 0:5:100
        , ylabel = "Percent"
        , title = column_names[column]
        )
        ylims!(ax, 0, 100)

        #Make violin plot with combined sample data
        violin!(fig[1,column]
        , repeat([x[column]]
        , size(metrics_file)[1]
        )
        , y[column]
        , show_median=true
        )
        next!(progress_bar_update)

        #Save file once all columns have been added
        if column == size(metrics_file)[2]
            save("Violin_metrics.png", fig)
        end
        next!(progress_bar_update)
    end
end

#Create violin plot of each sample's miRNA counts
function ViolinPlotMiRNA(miRNA_Counts_Files :: Vector{String})

    number_of_records :: Int64 = 4
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Creating miRNA Counts Violin Plots...")

    #Empty dictionary which will hold miRNA counts info from all samples
    miRNA_counts_dictionary = Dict{String, Vector{Int64}}()
    sample_names = Vector{String}()

    #Create dictionary with the sample name and miRNA counts as a key:value pair
    for file in miRNA_Counts_Files
        miRNA_file :: DataFrame = CSV.read(file, DataFrame)
        sample_name = first(split(file, "_"))
        counts = miRNA_file[!, :count]
        miRNA_counts_dictionary[sample_name] = counts
        append!(sample_names, [sample_name])
    end
    
    #Empty lists which will contain values for x and y axes
    #x-axis is vector of integers from 1 to the number of samples
    #y-axis is vector of vectors containing the miRNA counts data from each sample
    x_values :: Vector{Int64} = []
    miRNA_values :: Vector{Int64} = []

    #Function to add each sample's counts data to a combined array for ploting
    miRNA_counts_column(key) = [append!(miRNA_values, data) for data in miRNA_counts_dictionary[key]]

    #Function creates a vector containing the sample number repeated n times
    #where n is the number of different miRNA present in each sample
    repeat_numbers(num_repeats, iteration) = repeat([iteration+=1], num_repeats)

    #Make empty figure
    fig = Figure(resolution = (1800, 1200))
    next!(progress_bar_update)

    #Array comprehensions used to setup the arrays containing the x and y-axis data
    y = [miRNA_counts_column(sample) for sample in sample_names]
    x = [append!(x_values, repeat_numbers(length(y[i]), j)) for (i,j) in enumerate(0:length(miRNA_counts_dictionary)-1)]
    
    #Create violin plot with data from all samples analyzed
    ax = Axis(fig[1, 1]
    , xticks = (1:length(sample_names)
    , sample_names)
    , yticks = 0:20:length(sort(y, by=length, rev=true)[1])
    , ylabel = "number of miRNA"
    , title = "miRNA counts"
    )
    ylims!(ax, 0, length(sort(y, by=length, rev=true)[1]) + 20)
    next!(progress_bar_update)

    #Make violin plot with combined sample data
    violin!(fig[1,1]
    , x_values
    , miRNA_values
    , show_median=true
    )
    next!(progress_bar_update)

    #Save file
    save("Violin_miRNA_Counts.png", fig)
    next!(progress_bar_update)
end

#Create a file with a list of unique species names given a miRBase mature or hairpin fasta
function GetGenusNames(miRBase_Fasta::String)
    miRNA_reference :: IOStream = open(miRBase_Fasta, "r")
    genus_dictionary = Dict{String, String}()

    for line in eachline(miRNA_reference)
        genus_abbreviation = match(r">\K\w+", line)
        genus_name = match(r">.+MI\w+\s\K\w+", line)
        if isnothing(genus_name) === false
            genus_dictionary[genus_name.match] = genus_abbreviation.match
        end
    end

    close(miRNA_reference)

    output_genus_file :: IOStream = open("mirBase_genus_list.txt", "w")
    for (genus, abbreviation) in genus_dictionary
        write(output_genus_file, genus * "," * abbreviation * "\n")
    end

    close(output_genus_file)
end

#=
Function takes in mirBase mature miRNA fasta file with data from all available
organisms and pulls out only the miRNA data pertaining to the target organism.
=#
function CreateTargetOrganismFastaFile(Organism_Name :: String)

    #Create file with all genuses with miRBase annotations
    GetGenusNames("mature.fa")

    #Sets the target or organism's genus and species for labeling the output file
    organism_genus_name :: String = first(split(Organism_Name, " "))
    organism_species_name :: String = last(split(Organism_Name, " "))
    target_genus_abbreviation :: String = ""

    #Opens the IO streams and sets the output file names for the fasta and bowtie2 index
    genus_file :: IOStream = open("mirBase_genus_list.txt", "r")
    input_fasta_file::IOStream = open("mature.fa", "r")
    output_fasta_file::IOStream = open(organism_genus_name * "_" * organism_species_name * ".fa", "w")
    output_fasta_file_name :: String = organism_genus_name * "_" * organism_species_name * ".fa"
    bowtie2_reference_name :: String = organism_genus_name * "_" * organism_species_name

    #Loops through file containing the genus of all the organisms with miRNA data
    #If the target organism is in that list, the genus abbreviation is stored for use below
    for line in eachline(genus_file)
        if occursin(lowercase(organism_genus_name), lowercase(line))
            target_genus_abbreviation = last(split(line, ","))
            break
        end
    end

    #Loops through mirBase mature miRNA fasta looking for the header and sequence information
    #for the target genus. If there's a match, that information is added to a new fasta file.
    for line in eachline(input_fasta_file)
        input_genus_abbreviation = match(r">\K\w+", line)
        if isnothing(input_genus_abbreviation) === false
            if target_genus_abbreviation == input_genus_abbreviation.match
                write(output_fasta_file, line, "\n")
                write(output_fasta_file, readline(input_fasta_file), "\n")
            end
        end
    end

    close(input_fasta_file)
    close(output_fasta_file)
    close(genus_file)

    #Uses the newly created single organism fasta to build a bowtie2 index; will be used
    #for alignment with mirUtils to determine miRNA count information.
    wait(run(pipeline(
    `bowtie2-build 
    $output_fasta_file_name 
    $bowtie2_reference_name`
    , devnull)
    , wait=false
    ))

    TrashRemoval("mirBase_genus_list.txt")

    return bowtie2_reference_name, target_genus_abbreviation
end

#=
Takes in an array of adapter trimmed fastq files, the chosen reference name, and
whether or not a custom index was created. The correct bowtie2 reference is chosen
based on if a non-human reference was used. Then the trimmed fastq files are looped through
and aligned to the single organism bowtie2 reference. The bowtie2 output SAM is converted
to a BAM file and that is input into mirUtils to generate miRNA counts.
=#
function miRNADiscoveryCalculation(Trimmed_Fastq_Files :: Vector{String}, Organism_Name :: String, Reference_Choice :: String)
    number_of_records :: Int64 = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, 2.5, "Calculating the number of miRNA present...")

    #Set bowtie2 reference depending on whether or not a non-human reference was used
    if Reference_Choice === "y"
        bowtie2_reference_name, target_genus_abbreviation  = CreateTargetOrganismFastaFile(Organism_Name)
    else
        bowtie2_reference_name, target_genus_abbreviation  = CreateTargetOrganismFastaFile(Organism_Name)
    end

    #Loops through each sample and calulates how many miRNA were captured
    for fastq_file in Trimmed_Fastq_Files
        sample_name :: SubString{String} = match(r"sub_\K\w+(?=\.)", fastq_file).match
        wait(run(pipeline(
        `bowtie2
        --threads 12
        -x $bowtie2_reference_name
        -U $fastq_file
        -S $sample_name.miRNA.sam`
        , devnull)
        ,wait=false
        ))

        run(pipeline(
        `samtools
        view
        -@ 12
        -b
        -o $sample_name.miRNA.bam
        $sample_name.miRNA.sam`
        , devnull
        ))

        #mirUtils generates data for the specific miRNA were captured
        wait(run(pipeline(
        `mirUtils
        mbaseMirStats
        --version=v22
        --organism=$target_genus_abbreviation
        $sample_name.miRNA.bam`
        , devnull)
        , wait=false
        ))
                
        next!(progress_bar_update)
    end

    miRNA_files :: Vector{String} = CaptureTargetFiles("hairpin.hist")

    return miRNA_files
end

#Counts the number of unique miRNA groups each sample aligned to
function PlotMiRNACounts(miRNA_Files :: Vector{String})
    number_of_records :: Int64 = length(miRNA_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .25, "Counting the miRNA in each sample...")

    for file in miRNA_Files
        #Import mirUtils output file containing miRNA information
        miRNA_file :: DataFrame = CSV.read(file, DataFrame)

        #Checks column containing miRNA counts and adds miRNA above a set threshold
        threshold_of_one = sum(miRNA_file[!,:count] .>= 1)
        threshold_of_three = sum(miRNA_file[!,:count] .>= 3)
        threshold_of_five = sum(miRNA_file[!,:count] .>= 5)
        threshold_of_ten = sum(miRNA_file[!,:count] .>= 10)

        #Captures the name and read count of each miRNA with at least 1 read
        miRNA_at_threshold_of_one = miRNA_file[miRNA_file[!,:count] .>= 1, [:name, :count]]

        sample_name :: String = first(split(file, "."))

        #Colors for the bars in the plot
        colors = [:grey88, :skyblue2, :peachpuff, :lightsalmon]

        #Set x and y-axis values; y-axis is a vector containing each sample's counted miRNAs 
        x, y = 1:4, [threshold_of_one, threshold_of_three, threshold_of_five, threshold_of_ten]

        fig = Figure(resolution = (1800, 1200))
        ax = Axis(fig[1, 1]
        , xticks = (1:4
        , ["Threshold 1", "Threshold 3", "Threshold 5", "Threshold 10"])
        , yticks = 0:20:threshold_of_one
        , ylabel = "Number of miRNA"
        , title = sample_name * " miRNA Counts"
        )
        ylims!(ax, 0, threshold_of_one + 5)

        #Make sample barplot
        barplot!(ax
        , x
        , y
        , strokewidth = 1
        , color = colors
        , strokecolor = :black
        , bar_labels = :y
        )
        save(sample_name * "_miRNA_counts.pdf", fig)

        #Write output file containing miRNA names and read count data
        CSV.write(sample_name * "_miRNA_counts.csv", miRNA_at_threshold_of_one[!, [:name, :count]])
        next!(progress_bar_update)
    end

    pdf_miRNA_files :: Vector{String} = CaptureTargetFiles("_miRNA_counts.pdf")

    #Merge output PDFs into one file
    run(`gs 
    -q 
    -dBATCH 
    -dNOPAUSE 
    -sDEVICE=pdfwrite 
    -sOutputFile=ComboSeq_miRNA_plots.pdf 
    $pdf_miRNA_files`
    )
    #Removes individual sample plots
    TrashRemoval(pdf_miRNA_files)

    miRNA_counts_files :: Vector{String} = CaptureTargetFiles("_miRNA_counts.csv")

    return miRNA_counts_files

end

#Counts the number of unique mRNA groups each sample aligned to
function PlotMRNACounts(mRNA_Files :: Vector{String})
    number_of_records :: Int64 = length(mRNA_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records, .25, "Counting the mRNA in each sample...")

    for file in mRNA_Files
        #Import mirUtils output file containing mRNA information
        mRNA_file :: DataFrame = CSV.read(file, DataFrame)

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

        sample_name :: SubString{String} = match(r"sub_\K\w+(?=\.)", file).match

        #Colors for the bars in the plot
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
        ax = Axis(fig[1, 1]
        , xticks = (1:5
        , ["Threshold 1", "Threshold 10", "Threshold 100", "Threshold 500", "Threshold 1000"])
        , yticks = 0:1000:threshold_of_one
        , ylabel = "Number of mRNA"
        , title = sample_name * " mRNA Counts"
        )
        ylims!(ax, 0, threshold_of_one + 500)

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

    pdf_mRNA_files :: Vector{String} = CaptureTargetFiles("_mRNA_counts.pdf")

    #Merge output PDFs into one file
    run(`gs 
    -q 
    -dBATCH 
    -dNOPAUSE 
    -sDEVICE=pdfwrite 
    -sOutputFile=ComboSeq_mRNA_plots.pdf 
    $pdf_mRNA_files`
    )
    #Removes individual sample plots
    TrashRemoval(pdf_mRNA_files)

end

#Spot remove unecessary files
function TrashRemoval(Files_to_Delete :: Vector{String})
    for file in Files_to_Delete
        rm(file)
    end
end

function TrashRemoval(File_to_Delete :: String)
    rm(File_to_Delete)
end

#Remove intermediate files
function GarbageCollection()
    Files_to_Delete = (
    Subsampled_Files = CaptureTargetFiles("sub_")
    ,Deduplication_Files = CaptureTargetFiles("deduplication_statistics")
    ,Bowtie2_Files = CaptureTargetFiles(".bt2")
    ,mirUtils_outputs_pt1 = CaptureTargetFiles(".miRNA.")
    ,mirUtils_outputs_pt2 = CaptureTargetFiles("v22_cluster")
    )

    for set_of_files in Files_to_Delete
        for file in set_of_files
            rm(file)
        end
    end

    json_files = CaptureTargetFiles(".json")
    TrashRemoval(json_files)

end

function main()
    Salmon_Reference :: String, Organism_Name, Reference_Choice :: String = CreateReferenceFile(need_reference)
    Fastqs :: Vector{String} = CaptureTargetFiles(".fastq.gz")
    Downsampled_Fastqs = DownsampleRawReads(Fastqs)
    Trimmed_Fastqs :: Vector{String} = TrimAdapters(Downsampled_Fastqs)
    Salmon_SAM_Files :: Vector{String}, Salmon_mRNA_Counts :: Vector{String} = AlignWithSalmon(Trimmed_Fastqs,Salmon_Reference)
    Salmon_BAM_Files :: Vector{String} = ConvertSAMToBAM(Salmon_SAM_Files)
    PlotMRNACounts(Salmon_mRNA_Counts)
    Length_Files :: Vector{String} = CalculateReadLengthDistribution(Salmon_BAM_Files)
    PlotFragmentLengths(Length_Files)
    miRNA_Files :: Vector{String} = miRNADiscoveryCalculation(Trimmed_Fastqs, Organism_Name, Reference_Choice)
    miRNA_Counts_Files :: Vector{String} = PlotMiRNACounts(miRNA_Files)
    ViolinPlotMiRNA(miRNA_Counts_Files)
    Metrics_File :: String = CalculateSalmonMetrics(Trimmed_Fastqs)
    PlotMetrics(Metrics_File)
    ViolinPlotMetrics(Metrics_File)
    GarbageCollection()
end

main()
