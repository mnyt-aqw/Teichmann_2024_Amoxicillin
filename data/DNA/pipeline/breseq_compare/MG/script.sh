gdtools COMPARE \
    --format TSV \
    --output compare_MG.tsv \
     --reference /storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000005845.2/genomic.gbff \
     *.gd

gdtools COMPARE \
    --format HTML \
    --output compare_MG.html \
     --reference /storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000005845.2/genomic.gbff \
     *.gd