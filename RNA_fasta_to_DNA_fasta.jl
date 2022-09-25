# Script converts an input RNA fasta file into one with DNA sequences
# Use as: julia RNA_fasta_to_DNA_fasta.jl <RNA fasta file>

function RNAToDNAMappings(base::Char)
    base == 'C' && return 'C'
    base == 'G' && return 'G'
    base == 'U' && return 'T'
    base == 'A' && return 'A'
    base == 'N' && return 'N'
    base == 'R' && return 'R'
    base == 'Y' && return 'Y'
    base == 'S' && return 'S'
    base == 'W' && return 'W'
    base == 'K' && return 'K'
    base == 'M' && return 'M'
    base == 'H' && return 'H'
    base == 'B' && return 'B'
    base == 'V' && return 'V'
    base == 'D' && return 'D'
    throw(ErrorException("Unknown nucleotide: $base"))
end

function ParseFastaFile(Input_Fasta::String)
    input_fasta_file::IOStream = open(Input_Fasta, "r")
    output_fasta_file::IOStream = open(Input_Fasta * "_DNA", "w")
    counter::Int64 = 1
    for line in eachline(input_fasta_file)
        first_character = match(r"^.{1}", line).match
        if first_character == ">"
            write(output_fasta_file, line, "\n")
        else
            DNA_sequence = map(RNAToDNAMappings::Function, line)
            write(output_fasta_file, DNA_sequence, "\n")
        end
        counter +=1
    end
    close(input_fasta_file)
    close(output_fasta_file)
    mv(Input_Fasta * "_DNA", Input_Fasta, force=true)
end

function main()
    ParseFastaFile(ARGS[1])
end

main()