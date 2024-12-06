#!/usr/bin/env nextflow
/*
############################################

Pipeline for Breseq data analysis

############################################
*/
nextflow.enable.dsl=2


ref_BW = "${PWD}/../../data/RNA/references/ncbi_dataset/data/GCF_000750555.1/genomic.gbff"
ref_MG = "${PWD}/../../data/RNA/references/ncbi_dataset/data/GCF_000005845.2/genomic.gbff"
params.genome_mg1655 = "${PWD}/../../data/RNA/references/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna" 
params.input_reads = "${PWD}/../../data/DNA/data_raw/*/S*R{1,2}.fq.gz"

list_format = ["html", "tsv"]

params.directory_out = "${PWD}/../../data/DNA/pipeline"

workflow {
    Channel
        .fromFilePairs(params.input_reads)
        .map { sample_id, files ->
            def ref_file = sample_id.contains("dinB") || sample_id.contains("katE") ? ref_BW : ref_MG
            tuple(sample_id, files, file(ref_file))
        }
        .set { reads_ch }

    TrimGalore(reads_ch)

    Breseq(TrimGalore.out.trimmed)

    MultiQC(TrimGalore.out.fastqc_files.collect())

    e_coli_index()
    copy_number(TrimGalore.out.trimmed, e_coli_index.out.index)
    copy_number_bootstrap(copy_number.out.hits.collect(), ref_MG)
    breseq_split()
    Breseq_COMPARE_MG(breseq_split.out.MG)
}

process TrimGalore{
    publishDir = "${params.directory_out}/TrimGalore/"
    cpus 10

    input:
    tuple val(sample_id), path(reads), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}*_R{1,2}_val_{1,2}.fq.gz"), path(ref), emit: trimmed
    path "*fastqc.zip", emit: fastqc_files

    script:
    """
    trim_galore --paired --fastqc --adapter "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" --phred33 -e 0.1 --quality 28 --cores ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

process Breseq {
    publishDir = "${params.directory_out}Breseq/"
    cpus 10

    input:
    tuple val(sample_id), path(reads), path(ref)

    output:
    path "${sample_id}/", emit: mutation_analysis

    script:
    """
    breseq -r ${ref} ${reads[0]} ${reads[1]} --num-processors ${task.cpus} --name ${sample_id} --polymorphism-prediction
    mkdir ${sample_id}
    mv 0* ${sample_id}/.
    mv data/ ${sample_id}/.
    mv output/ ${sample_id}/.
    """
}

process MultiQC {
    publishDir = "${params.directory_out}/MultiQC/"
    cpus 1

	input:
	path files

	output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc -f *fastqc.zip
    """
}

process ReplaceGenomeID_mg1655 {
    input:
    tuple val(id) path(file)

    output:
    path "modified_${sampleId}.gd", emit: modified

    script:
    """
    sed 's/NC_000913/NC_000913.3/g' ${file} > modified_S${id}.gd
    """
}

process e_coli_index {
    cpus 5

	output:
    path "MG1655*", emit: index

    script:
    """
    bowtie2-build ${params.genome_mg1655} MG1655
    """
}

process copy_number {
    publishDir = "${params.directory_out}/Copy_number/"
    cpus 5

	input:
    tuple val(sample_id), path(reads), path(gbff)
    path index

	output:
    path "${sample_id}.hits", emit: hits
    path "${sample_id}_align_summary.txt"

    script:
    """
    bowtie2 -p ${task.cpus}  -x MG1655 -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam 2> ${sample_id}_align_summary.txt
    samtools view -@ ${task.cpus} -Sb  ${sample_id}.sam > ${sample_id}.bam
    samtools rmdup -S ${sample_id}.bam  ${sample_id}_rmdup.bam
    samtools view ${sample_id}_rmdup.bam | perl -lane 'print "\$F[2]\\t\$F[3]"' > ${sample_id}.hits
    sed -i "s/NC_000913\\.3/NC_000913/g" ${sample_id}.hits
    sed -i "/\\*/d" ${sample_id}.hits
    rm ${sample_id}.sam
    """
}

process copy_number_bootstrap {
    publishDir  "${params.directory_out}/Copy_number/", mode: "copy"
    cpus 1
    conda = "/storage/marwe/envs/lisa_RNA/"

    input:
    path hit_files 
    path gbkfile

    output:
    path "gene_cp_nr.csv"

    script:
    """
    Rscript --vanilla ${PWD}/nextflow_scripts/cp_nr.r
    """
}

process breseq_split {
    publishDir  "${params.directory_out}/Copy_number/", mode: "copy"
    cpus 1

    output:
    path("MG/*.gd"), emit: MG
    path("BW/*.gd"), emit: BW

    script:
    """
    #!/usr/bin/env python
    import os
    import glob
    import shutil
    import subprocess

    # Set the main directory path
    main_dir = "/storage/marwe/LisaRNA_DNA/data/DNA/pipeline/breseq_compare/Breseq/"
    our_dir = "./"

    # Set the path to the reference GenBank file
    ref_mg_gbk = "/storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000005845.2/genomic.gbff"
    ref_bw_gbk = "/storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000750555.1/genomic.gbff"

    # Create the output directories if they don't exist
    os.makedirs(os.path.join(our_dir, "MG"), exist_ok=True)
    os.makedirs(os.path.join(our_dir, "BW"), exist_ok=True)

    for root, dirs, files in os.walk(main_dir):
        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            output_gd_path = os.path.join(dir_path, "data", "annotated.gd")

            if "MG1655" in output_gd_path:
                destination = os.path.join(our_dir, "MG", f"{dir_name}.gd")
                shutil.copyfile(output_gd_path, destination)
            else:
                destination = os.path.join(our_dir, "BW", f"{dir_name}.gd")
                shutil.copyfile(output_gd_path, destination)

    """
}