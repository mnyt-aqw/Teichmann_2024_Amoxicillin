# RNA analysis
___

## Pre processing

The Nextflow pipeline is for RNA-seq data analysis and consists of several processes. It starts by reading in the raw read files through the `read_ch` channel. The `TrimGalore` process trims the read files and generates quality control outputs, including trimmed reads and quality reports, which are then processed by the `MultiQC` process to produce a single quality control report.

The reference genome is then downloaded using the `Download_ref` process, and a `Bowtie2` index is created using the `Bowtie_index` process for aligning the trimmed reads to the reference genome. The `Bowtie_Samtools` process aligns the trimmed reads to the reference genome and converts the output to a sorted binary format (BAM). The resulting BAM files are processed by the `FeatureCounts` process to generate read count files, which are then processed by the `Add_gene_names` process to add gene names to the read counts.

Each process has a defined input and output, along with a script that performs the respective step in the analysis. For example, the `TrimGalore` process takes as input the file name and path of the raw reads and generates trimmed reads, quality control outputs, and trimming reports. The `FeatureCounts` process takes as input the BAM files and reference genome and outputs read count files. The output of each process is saved in a specified directory using the `publishDir` command.

Note that some of the processes, like `Bowtie_Samtools` and `Add_gene_names` reference a Python script in their script commands. These scripts contain additional code and details for the respective steps in the analysis, beyond what is specified in the Nextflow pipeline.

This Nextflow pipeline is for RNA-seq data analysis and consists of the following processes:

## Pipeline steps 

1. Trimming: The raw RNA-seq reads are trimmed using TrimGalore to remove low quality bases and adapters.
2. Quality control: The quality of the trimmed reads is assessed using MultiQC, which generates a single quality control report.
3. Reference download: The reference genome and annotations are downloaded from NCBI and unzipped.
4. Indexing: The reference genome is indexed using Bowtie2, which prepares it for mapping.
5. Mapping: The trimmed reads are mapped to the reference genome using Bowtie_Samtools.
6. Quantification: Gene expression levels are quantified using FeatureCounts, which counts the number of reads that map to each gene.
7. Adding gene names: The gene expression quantification results are annotated with gene names using Add_gene_names.
