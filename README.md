# Combo-Seq-Pipeline
Basic analysis of NEXTFLEX Combo-Seq libraries

[![CI](https://github.com/cookienocreams/combo_seq/actions/workflows/CI.yaml/badge.svg)](https://github.com/cookienocreams/combo_seq/actions/workflows/CI.yaml)

## Script dependencies
- bowtie2
- samtools
- cutadapt
- salmon

 ## Note: This currently only compiles under the Upcoming release: v1.10.0-beta2
 There is a bug in the Plots library that prevents packages from being complied using PackageCompiler, see [Plots.lj issue #825](https://github.com/JuliaLang/PackageCompiler.jl/issues/825).
 It should be fixed in v1.10, once that releases. Versions < v1.8.5 may also work.
 
## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/Combo-Seq-Pipeline.git combo_seq
cd combo_seq
```
Activate the downloaded `julia` environment.
```julia
using Pkg
Pkg.activate("./")
```
The next step is to install all libraries and their dependencies.
```julia
Pkg.instantiate()
```

The last step is to create the precompiled executable. Make sure to set the correct paths for your machine.

```julia
using PackageCompiler
PackageCompiler.create_app("./", "/home/user/combo_seq_app", incremental=true, precompile_execution_file="./src/combo_seq.jl", include_lazy_artifacts=true)
```

There are a few options that can be changed if desired. Use `-h` or `--help` to see help with descriptions of options.

```bash
/home/user/combo_seq_app/bin/combo_seq --help
```

The app can be run using the `combo_seq` executable in the `/home/user/combo_seq_app/bin` folder in a folder containing fastqs to be analyzed.

## Basic usage

In a folder containing Combo-Seq fastq files and a bowtie2 miRNA reference `/home/user/miRNA` and a Salmon reference `/home/user/hg38_mRNA` use:

```julia
/home/user/combo_seq_app/bin/combo_seq --mrna /home/user/hg38_mRNA --mirna /home/user/miRNA
```

If you do not have a Salmon reference already made and would like the program to make one, you can pass a link to both a transcriptome and genome fasta file. 
Note this can take a while for large genomes.

```julia
/home/user/combo_seq_app/bin/combo_seq --need-reference \
                                       --transcript https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz \
                                       --genome https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz
```

That will make a Salmon reference called `hsa` located in the `data` folder where the program was run. It can be used in subsequent runs as follows:

```julia
/home/user/combo_seq_app/bin/combo_seq --mrna /home/user/combo_seq/data/hsa
```
