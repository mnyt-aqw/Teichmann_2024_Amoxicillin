gdtools COMPARE \
    --format TSV \
    --output compare_BW.tsv \
     --reference /storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000750555.1/genomic.gbff \
     *.gd

gdtools COMPARE \
    --format HTML \
    --output compare_BW.html \
     --reference /storage/marwe/LisaRNA_DNA/data/RNA/references/ncbi_dataset/data/GCF_000750555.1/genomic.gbff \
     *.gd