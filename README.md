# Differences between antibiotics affecting the resistance evolution trajectories of Escherichia coli

The `scripts` directory contains scripts for running processes. The `documentation`contains workflow notebooks and other information about how the scripts were run.

## DNA

The pipeline for the mutation analysis and copy number can be found here: `scripts/nextflow_RNA/main.nf`.

The code for generating the Breseq tsv and html files can be found here `data/DNA/pipeline/breseq_compare/*/script.sh`

## RNA

The transcriptomics pipeline can be found here: `scripts/nextflow_RNA/main.nf`

## Data analysis

The notebook with the scientific questions can be found here: `documentation/RNA/scientific_questions/Scientific_questions_notebook.ipynb`.

## Container

This directory contains 3 relevant files.

`pull_container.sh`: Will pull the data analysis container from DockerHub.

`launch_jupyter.sh`: This will launch a jupyter server from the container. Meaning that all relevant data analysis packages are there in the right versions.

`render_questions.sh`: This will render the `Scientific_questions_notebook.ipynb`notebook into an interactive html file.

## Method flow chart

This is the bioinformatics workflow.  

```mermaid
graph TD
    subgraph Transcriptomics
        A1[Wetlab work] --> RNA[RNA-seq]
        RNA --> B1[TrimGalore: Trimming & QC]
        B1 --> B2[MultiQC: Quality Control]
        B1 --> B3[Download MG1655 References]
        B3 --> B4[Bowtie2: Index Building]
        B1 --> B5[Bowtie2 & Samtools: Alignment]
        B4 --> B5
        B5 --> B6[FeatureCounts: Read Quantification]
    end

    subgraph WGS Analysis
        A2[Wetlab work] --> DNA[DNA-seq]
        DNA --> C1[TrimGalore: Trimming & QC]
        C1 --> C2[MultiQC: Quality Control]
        C1 --> C3[Breseq: mutation screen]
        C4[E. coli Index: MG and BW]
        C4 --> C5[Bowtie2 & CNOGpro: Copy Number Analysis]
        C1 --> C5
        C4 --> C3
    end

    subgraph Jupyter Notebook Analysis
        B6 --> D2[edgeR: Differential Expression Analysis]
        D2 --> D4[GO: Gene Ontology Analysis]
        D2 --> D5[KEGG: Pathway Analysis]
        C3 --> D6[Mutation Analysis]
        C5 --> D7[Copy Number Analysis]
        D2 --> D8[Gene: transcript analysis]
    end
```

## Which environment was used?  

All data processing done using Nextflow pipelines or scripts was done using a Conda environment which you can find here `container/conda_env.yml`. The Jupyter notebook data analysis was performed using the `container/ecoli_evolution.sif` container which can be pulled using the `container/pull_container.sh` script. You can also build it using the `container/Dockerfile` script.