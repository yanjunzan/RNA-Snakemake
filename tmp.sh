STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --sjdbGTFtagExonParentTranscript Parent \
    --genomeDir /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/ \
    --genomeFastaFiles /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/Barley_MorexV3_pseudomolecules.fasta \
    --sjdbGTFfile  /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.gff3\
    --sjdbOverhang 99 \
    --alignIntronMax 20000
    hisat2 -p 2 --met-file /Users/yanjunzan/Documents/projects/barley_en/results/logs/ERR4393948_met.txt -x /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/Barley_MorexV3_pseudomolecules  -1 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/ERR4393948_R1_trimmed.fq.gz -2 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/ERR4393948_R2_trimmed.fq.gz | samtools view -Sb -F 4 -o /Users/yanjunzan/Documents/projects/barley_en/results/mapped/ERR4393948.bam
    /Users/yanjunzan/OneDrive - Sveriges Lantbruksuniversitet/Projects/barley_en/bin/RNA_seq/Snakefile:
    hisat2 -p 1 --met-file /Users/yanjunzan/Documents/projects/barley_en/results/logs/SRR8992621_met.txt -x /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/Barley_MorexV3_pseudomolecules             -1 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/SRR8992621_R1_trimmed.fq.gz -2 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/SRR8992621_R2_trimmed.fq.gz > /Users/yanjunzan/Documents/projects/barley_en/results/mapped/SRR8992621.sam
