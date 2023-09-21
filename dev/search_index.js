var documenterSearchIndex = {"docs":
[{"location":"#NEXTFLEX-Combo-Seq","page":"NEXTFLEX Combo-Seq","title":"NEXTFLEX Combo-Seq","text":"","category":"section"},{"location":"#Index","page":"NEXTFLEX Combo-Seq","title":"Index","text":"","category":"section"},{"location":"","page":"NEXTFLEX Combo-Seq","title":"NEXTFLEX Combo-Seq","text":"","category":"page"},{"location":"#Functions","page":"NEXTFLEX Combo-Seq","title":"Functions","text":"","category":"section"},{"location":"","page":"NEXTFLEX Combo-Seq","title":"NEXTFLEX Combo-Seq","text":"Documentation for combo_seq.jl","category":"page"},{"location":"","page":"NEXTFLEX Combo-Seq","title":"NEXTFLEX Combo-Seq","text":"MissingReferenceError\r\nIncorrectOrganismError\r\nMissingFastqsError\r\nConfig\r\nget_genus_names\r\ncreate_target_organism_fasta_file\r\ncapture_target_files\r\nprogress_bar_update\r\nget_read_q_score!\r\nmake_bowtie2_reference\r\nmake_salmon_reference\r\ncreate_reference_file\r\nparse_fastq_files\r\ntrim_adapters\r\nalign_with_salmon\r\ncalculate_salmon_metrics\r\ngenerate_mirna_counts\r\nmake_metrics_violin_plot\r\nmirna_discovery_calculation\r\nplot_mirna_counts\r\ncalculate_read_length_distribution\r\nplot_fragment_lengths\r\nplot_metrics\r\nfind_common_rnas\r\nwrite_common_rna_file\r\nplot_clustering\r\nremove_files\r\nremove_intermediate_files","category":"page"},{"location":"#combo_seq.MissingReferenceError","page":"NEXTFLEX Combo-Seq","title":"combo_seq.MissingReferenceError","text":"struct MissingReferenceError\n    message::String\nend\n\nA structure that represents a custom exception for a missing reference.\n\n\n\n\n\n","category":"type"},{"location":"#combo_seq.IncorrectOrganismError","page":"NEXTFLEX Combo-Seq","title":"combo_seq.IncorrectOrganismError","text":"struct IncorrectOrganismError\n    message::String\nend\n\nA structure that represents a custom exception for an incorrectly specified orgamism name.\n\n\n\n\n\n","category":"type"},{"location":"#combo_seq.MissingFastqsError","page":"NEXTFLEX Combo-Seq","title":"combo_seq.MissingFastqsError","text":"struct MissingFastqsError\n    message::String\nend\n\nA structure that represents a custom exception for missing fastq files.\n\n\n\n\n\n","category":"type"},{"location":"#combo_seq.Config","page":"NEXTFLEX Combo-Seq","title":"combo_seq.Config","text":"struct Config\n    need_reference::Bool\n    mrna::Union{String, Nothing}\n    fasta::Union{String, Nothing}\n    transcript::Union{String, Nothing}\n    genome::Union{String, Nothing}\n    organism::Union{String, Nothing}\n    threads::Int\nend\n\nA structure that stores the command-line arguments.\n\nFunctionality\n\nContains a check to ensure that a Salmon and miRNA fasta file are passed as arguments. \n\nIf they were not, an error will be returned. This is to reduce the chances of an alignment  errors downstream.\n\n\n\n\n\n","category":"type"},{"location":"#combo_seq.get_genus_names","page":"NEXTFLEX Combo-Seq","title":"combo_seq.get_genus_names","text":"Create a file with a list of unique species names given a miRNA mature or hairpin fasta.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.create_target_organism_fasta_file","page":"NEXTFLEX Combo-Seq","title":"combo_seq.create_target_organism_fasta_file","text":"Function takes in mirBase mature miRNA fasta file with data from all available organisms and pulls out only the miRNA data pertaining to the target organism.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.capture_target_files","page":"NEXTFLEX Combo-Seq","title":"combo_seq.capture_target_files","text":"capture_target_files(files_to_capture::AbstractString, directory::AbstractString=\".\")\n\nList all files in a directory. \n\nCheck to see if each file contains the target file's string. \n\nExample\n\njulia> capture_target_files(\".txt\")\n3-element Vector{String}:\n \"file1.txt\"\n \"file2.txt\"\n \"file3.txt\"\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.progress_bar_update","page":"NEXTFLEX Combo-Seq","title":"combo_seq.progress_bar_update","text":"Update progress bar on the command line.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.get_read_q_score!","page":"NEXTFLEX Combo-Seq","title":"combo_seq.get_read_q_score!","text":"get_read_q_score!(line::String, q_score_list::Vector{Number})\n\nCalculate the quality score for a given quality string.\n\nThis function accepts ASCII-encoded quality score and produces the average q score  using a scalar value that estimates the 'total error probability' for that record. It means that if you want to calculate average quality score, you don't just sum up all the phred scores and find the average, but rather consider the  probability distribution of the scores.\n\nThe Phred score for a base Q is calculated as Q = ASCII value of quality score - 33. The error probability P for that base is then calculated as P = 10^(-Q/10). Then these probabilities are summed up for all the bases to derive total error.\n\nArguments\n\nq_score - The quality scores of the current line. q_score_list - The vector containing all calculated quality scores.\n\nReturns\n\nThe updated q score list the line's average quality score.\n\nReference\n\nIllumina's explanation on quality scores and their corresponding error rates: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html\n\nExample\n\njulia> get_read_q_score!(\"FGFEEE<FC<FGGGGGGGGGGGFFGFG8<@8CFFGF8EFGGCGFGGGGGGFG\", [36.2, 35.9])\n3-element Vector{Float64}:\n 36.2\n 35.9\n 35.7\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.make_bowtie2_reference","page":"NEXTFLEX Combo-Seq","title":"combo_seq.make_bowtie2_reference","text":"Make miRNA bowtie2 reference to align fastqs with.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.make_salmon_reference","page":"NEXTFLEX Combo-Seq","title":"combo_seq.make_salmon_reference","text":"Make Salmon reference for mRNA alignment.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.create_reference_file","page":"NEXTFLEX Combo-Seq","title":"combo_seq.create_reference_file","text":"Function asks the user if they would like to download a non-human reference for use by Salmon during the alignment step. If yes, the files are downloaded and the genome names are extracted to name the Salmon directory.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.parse_fastq_files","page":"NEXTFLEX Combo-Seq","title":"combo_seq.parse_fastq_files","text":"parse_fastq_files(fastqs::Vector{String}, sample_names::Vector{SubString{String}})\n\nCalculate read counts, percent dimers, and average Q scores for a given set of fastq files.\n\nParameters\n\nfastqs: A vector of filenames for the fastq files.\nsample_names: A vector of sample names corresponding to the fastq files.\n\nReturns\n\nThree dictionaries containing read counts, percent dimers, and average Q scores  indexed by sample name.\n\nExample\n\njulia> parse_fastq_files([\"sample1.fastq.gz\", \"sample2.fastq.gz\", \"sample3.fastq.gz\"], \n                        [\"sample1\", \"sample2\", \"sample3\"])\n# Returns:\n# - A dictionary of read counts.\n# - A dictionary of percent dimers.\n# - A dictionary of average Q scores.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.trim_adapters","page":"NEXTFLEX Combo-Seq","title":"combo_seq.trim_adapters","text":"trim_adapters(fastqs::Vector{String}, sample_names::Vector{SubString{String}})\n\nTrim the 3' adapter from each read in the provided fastq files.\n\nCombo-Seq Adapters:\n\n5' Adapter: GATCGTCGGACTGTAGAACTCTGAACNNNN \n\n3' Adapter: AAAAAAAAAA\n\n3' bases with a quality score below 20 are trimmed. Reads shorter than 16 bases or  those that weren't trimmed are discarded.\n\nParameters\n\nfastqs: A vector of filenames for the fastq files.\nsample_names: A vector of sample names corresponding to the fastq files.\n\nReturns\n\nA list of filenames for the trimmed fastq files.\n\nExample\n\njulia> trim_adapters([\"sample1.fastq.gz\", \"sample2.fastq.gz\", \"sample3.fastq.gz\"], \n                    [\"sample1\", \"sample2\", \"sample3\"])\n# Returns:\n[\"sample1.cut.fastq\", \"sample2.cut.fastq\", \"sample3.cut.fastq\"]\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.align_with_salmon","page":"NEXTFLEX Combo-Seq","title":"combo_seq.align_with_salmon","text":"align_with_salmon(trimmed_fastq_files::Vector{String}\n                    , salmon_reference::AbstractString\n                    , sample_names::Vector{SubString{String}}\n                    , need_reference::Bool\n                    )\n\nThis function aligns trimmed fASTQ files to a reference genome using Salmon.\n\nThe library type (libType) is set to stranded forward (SF) since the Combo-Seq protocol is stranded in that direction.\n\nArguments\n\ntrimmed_fastq_files: A vector of strings containing the paths to the trimmed fASTQ files.\nsalmon_reference: A string containing the path to the Salmon reference genome.\nsample_names: A vector of strings containing the names of the samples.\nneed_reference: A string indicating whether or not a custom Salmon index was created.\n\nReturns\n\nThe function returns a vector of strings containing the paths to the Salmon output files.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.calculate_salmon_metrics","page":"NEXTFLEX Combo-Seq","title":"combo_seq.calculate_salmon_metrics","text":"calculate_salmon_metrics(trimmed_fastq_files::Vector{String}\n                        , read_count_dict::Dict{String, Int64}\n                        , sample_names::Vector{SubString{String}}\n                        , dimer_count_dict::Dict{String, Float64}\n                        , q_score_dict::Dict{String, Number}\n                        )\n\nAlign fastq with specified reference RNA types.\n\nDivide reads aligned with total reads to calculate percent alignment. Store  alignment information in a dictionary.\n\n#Example\n\njulia> calculate_salmon_metrics(\"sample1.cut.fastq\"\n                        , Dict(\"sample1\" => 6390314, \"sample2\" => 5000000, \"sample3\" => 7052928))\n# Returns:\n# - \"ComboSeq_metrics.csv\"\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.generate_mirna_counts","page":"NEXTFLEX Combo-Seq","title":"combo_seq.generate_mirna_counts","text":"generate_mirna_counts(input_sam_file::AbstractString)\n\nThis function creates a DataFrame containing all aligned miRNA and the number of times  they appear.\n\nReturns\n\nA DataFrame with the following columns:\n\nname: The name of the miRNA.\ncount: The number of times the miRNA appears in the SAM file.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.make_metrics_violin_plot","page":"NEXTFLEX Combo-Seq","title":"combo_seq.make_metrics_violin_plot","text":"Create violin plot of each sample's RNA aligment metrics\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.mirna_discovery_calculation","page":"NEXTFLEX Combo-Seq","title":"combo_seq.mirna_discovery_calculation","text":"mirna_discovery_calculation(trimmed_fastq_files::Vector{String}\n                            , sample_names::Vector{SubString{String}}\n                            , bowtie2_reference_name::AbstractString\n                            , organism_name::AbstractString\n                            , need_reference::AbstractString\n                            )\n\nAlign trimmed fastq files to the single organism miRNA bowtie2 reference. \n\nThe bowtie2 output SAM is input into function to generate miRNA read counts.\n\n#Example\n\njulia> mirna_discovery_calculation([\"sample1.cut.fastq\",\"sample2.cut.fastq\",\"sample3.cut.fastq\"]\n                                , [\"sample1\", \"sample2\", \"sample3\"])\n3-element Vector{DataFrame}:\n469x2 DataFrame\n Row │ name             count\n     │ String           Int64\n─────┼────────────────────────\n   1 │ hsa-miR-1307-3p      1\n   2 │ hsa-miR-425-5p      58\n   3 │ hsa-miR-488-5p       2\n [...]\n\n 3-element Vector{String}:\n \"sample1.miRNA.sam\"\n \"sample2.miRNA.sam\"\n \"sample3.miRNA.sam\"\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.plot_mirna_counts","page":"NEXTFLEX Combo-Seq","title":"combo_seq.plot_mirna_counts","text":"plot_mirna_counts(mirna_counts_dfs::Vector{DataFrame}\n                    , sample_names::Vector{SubString{String}}\n                    )\n\nPlot counts of the number of unique miRNA each sample aligned to.\n\n#Example\n\njulia> plot_mirna_counts([sample1_counts_df,sample2_counts_df,sample3_counts_df]\n                        , [\"sample1\", \"sample2\", \"sample3\"])\n3-element Vector{String}:\n \"sample1_miRNA_counts.csv\"\n \"sample2_miRNA_counts.csv\"\n \"sample3_miRNA_counts.csv\"\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.calculate_read_length_distribution","page":"NEXTFLEX Combo-Seq","title":"combo_seq.calculate_read_length_distribution","text":"Use aligned SAM files to create a read length distribution file.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.plot_fragment_lengths","page":"NEXTFLEX Combo-Seq","title":"combo_seq.plot_fragment_lengths","text":"Create barplot of each sample's fragment lengths based on the number of reads found at each length\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.plot_metrics","page":"NEXTFLEX Combo-Seq","title":"combo_seq.plot_metrics","text":"Create barplot of each sample's RNA aligment metrics.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.find_common_rnas","page":"NEXTFLEX Combo-Seq","title":"combo_seq.find_common_rnas","text":"find_common_rnas(counts_files::Vector{String}\n                , sample_names::Vector{SubString{String}} \n                ; rna_type::String=\"miRNA\", read_count_dict::Dict{String\n                , Int64}=Dict{String, Int64}()\n                )\n\nReturn the RNAs present in all samples, together with their counts and reads per million reads (RPM) or transcripts per million (TPM).\n\nArguments\n\ncounts_files: A vector of files containing sample RNA counts.\nsample_names: A vector of sample names.\nrna_type: A string specifying the type of RNA (\"miRNA\" or \"mRNA\"). Default is \"miRNA\".\nread_count_dict: (Optional for \"mRNA\" type) A dictionary where keys are sample names and \n\nvalues are the total number of reads  for the corresponding sample. This is required for \"miRNA\" type to compute RPM values.\n\nReturns\n\nA Tuple of three items:\nA list of RNAs that are common to all samples.\nA list of dictionaries where each dictionary contains the RNA counts for one sample.\nA list of dictionaries where each dictionary contains the RPM (for miRNA) or TPM \n(for mRNA) values for the RNAs in one sample.\n\nExample\n\nFor miRNA:\n\njulia> find_common_rnas([\"sample1_miRNA_counts.csv\", \"sample2_miRNA_counts.csv\"],\n                        [\"sample1\", \"sample2\"],\n                        rna_type=\"miRNA\",\n                        read_count_dict=Dict(\"sample1\" => 6390314, \"sample2\" => 5000000))\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.write_common_rna_file","page":"NEXTFLEX Combo-Seq","title":"combo_seq.write_common_rna_file","text":"write_common_rna_file(rna_names::Base.KeySet{String, Dict{String, Int64}},\n                      rna_info::Vector{Dict{String, Float64}},\n                      normalization_info::Vector{Dict{String, Float64}},\n                      sample_names::Vector{SubString{String}};\n                      rna_type::String=\"miRNA\")\n\nWrite the common RNAs found in all samples to an output file. This function generates two sets of files for each sample:  one containing raw RNA counts and another containing normalized values (RPM for miRNA or TPM for mRNA). It then  combines these files for all samples to produce comprehensive outputs, and subsequently conducts a principal component  analysis and UMAP on the common RNA data.\n\nArguments\n\nrna_names: A set containing the names of RNAs that are common across all samples.\nrna_info: A vector of dictionaries with each dictionary containing the RNA counts for one sample.\nnormalization_info: A vector of dictionaries, each containing the RPM (for miRNA) or TPM (for mRNA) values for the RNAs in one sample.\nsample_names: A vector of sample names.\nrna_type: (Optional) A string specifying the type of RNA (\"miRNA\" or \"mRNA\"). Default is \"miRNA\".\n\nOutputs\n\nGenerates files for each sample named like sample_common_miRNAs.tsv and sample_common_miRNAs_RPM.tsv (or their mRNA equivalents).\nCombines these files to produce comprehensive outputs named Common_miRNAs.tsv and Common_RPM_miRNAs.tsv (or their mRNA equivalents).\nConducts principal component analysis and UMAP based on common RNA data.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.plot_clustering","page":"NEXTFLEX Combo-Seq","title":"combo_seq.plot_clustering","text":"plot_clustering(Common_miRNA_File::String)\n\nGenerates Principal Component Analysis (PCA) and Uniform Manifold Approximation and  Projection (UMAP) plots for libraries given a common RNA file.\n\nThis function performs the following steps:\n\nReads the common RNA file and creates a DataFrame.\nExtracts sample names and [mi/m]RNA names.\nEnsures that there are at least two samples and more than one miRNA present in the dataset.\nTransforms the data using PCA and prepares it for plotting.\nIf there are at least two principal components, plots the PCA and saves it as \"common[mi/m]RNAPCA.png\".\nIf there are at least three principal components, performs UMAP dimensionality reduction and clustering using K-medoids.\nPlots the UMAP and saves it as \"common[mi/m]RNAUMAP.png\".\nWrites the PCA and UMAP information to separate CSV files for cluster tracking.\n\nArguments\n\nCommon_[mi/m]RNA_File::String: Path to the common RNA RPM/TPM file.\n\nOutputs\n\nCreates and saves PCA and UMAP plots as \"common[mi/m]RNAPCA.png\" and \"common[mi/m]RNAUMAP.png\".\nWrites PCA and UMAP information to \"PCAinformation.csv\" and \"UMAPinformation.csv\".\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.remove_files","page":"NEXTFLEX Combo-Seq","title":"combo_seq.remove_files","text":"Spot remove unecessary intermediate files.\n\n\n\n\n\n","category":"function"},{"location":"#combo_seq.remove_intermediate_files","page":"NEXTFLEX Combo-Seq","title":"combo_seq.remove_intermediate_files","text":"Remove all intermediate files.\n\n\n\n\n\n","category":"function"}]
}
