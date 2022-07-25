# download, # fastqc # trim # hisat # index # sort # BAM qc #
# /Users/yanjunzan/anaconda3/envs/Stripes/etc/asperaweb_id_dsa.openssh
#/Users/yanjunzan/anaconda3/envs/Stripes/bin/aspera
#print(config["ref_dir"])
import os
import pandas as pd
configfile: "./config.yaml"

###########################
# define sample
###########################

df = pd.read_csv(config["sra_list"],sep="\t")
SAMPLES = df["Run"].tolist()


WORKING_DIR = config["result_dir"]
RESULT_DIR = config["result_dir"]

###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function checks if raeds are single or paired end"""
    if os.stat(sample).st_size < 100:
        return True
    else:
        return False



rule all:
    input:
        #fw  = expand(WORKING_DIR + "fastq/{sample}_1.fastq",sample = SAMPLES),
        qcfile = expand(RESULT_DIR + "fastp/{sample}.html",sample=SAMPLES),
        fq1 = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",sample = SAMPLES),
        fq2 = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",sample = SAMPLES),
        bam = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES)
    message:
        "Job done! Removing temporary directory"

#
# rule get_sra_paired:
#     conda:
#         "./mapping.yaml"
#     output:
#         expand("../../results/SRA_download/{sample_ids}.sra",sample_ids=samples)
rule get_SRR_files:
    output:
        fw = temp(WORKING_DIR + "fastq/{sample}_1.fastq"),
        rev= temp(WORKING_DIR + "fastq/{sample}_2.fastq")
    conda:
        "./mapping.yaml"
    params:
       SRA = "{sample}",
       DIR = config["result_dir"]+"fastq/"
    message:
        "using fastq-dump to download SRA data files to {output.fw}."
    conda:
        "./mapping.yaml"
    shell:
        "touch {output.rev}; fastq-dump --split-files {params.SRA} -O {params.DIR}"


rule fastp:
    input:
        fw = WORKING_DIR + "fastq/{sample}_1.fastq",
        rev= WORKING_DIR + "fastq/{sample}_2.fastq"
    output:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    conda:
        "./mapping.yaml"
    message:"trimming {wildcards.sample} reads to {output.fq1}"
    threads: 1
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        #paired = "{lambda wildcards:paired[wildcards.SAMPLES]}",
        #sampleName = "fastq/{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    # run:
    #     if sample_is_single_end(input.rev):
    #         shell("fastp --thread {threads}  --html {output.html} \
    #         --qualified_quality_phred {params.qualified_quality_phred} \
    #         --in1 {input.fw} --out1 {output} \
    #         2> {log}; \
    #         touch {output.fq2}")
    #     else:
    #         shell("fastp --thread {threads}  --html {output.html} \
    #         --qualified_quality_phred {params.qualified_quality_phred} \
    #         --in1 {input.fw} --in2 {input.rev} --out1 {output.fq1} --out2 {output.fq2}; \
    #         2> {log}")
    shell:
        "fastp --thread {threads}  --html {output.html} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        --in1 {input.fw} --in2 {input.rev} --out1 {output.fq1} --out2 {output.fq2}; \
        2> {log}"

rule index_hisat:
    conda:
        "mapping.yaml"
    input:
        config["working_dir"] +"MorexV3/"+ config["ref"]
    output:
        [config["working_dir"]+ "MorexV3/" + config["ref"].replace("fasta","") + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome {output}"
    params:
        index  = config["working_dir"] + "MorexV3/" + config["ref"].replace("fasta",""),
        genome = config["working_dir"] + "MorexV3/" + config["ref"]
    threads: 1
    shell:
        "hisat2-build -p {threads} {params.genome} {params.index} --quiet"

rule hisat_mapping:
    conda:
        "mapping.yaml"
    input:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        indexFiles = [config["working_dir"]+ "MorexV3/" + config["ref"].replace("fasta","") + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = protected(WORKING_DIR + "mapped/{sample}.bam"),
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = config["working_dir"]+ "MorexV3/" + config["ref"].replace(".fasta",""),
        sampleName = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files {params.sampleName}."
    threads: 2
    # run:
    #     if sample_is_single_end(params.sampleName):
    #         shell("hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
    #         -U {input.fq1} > {output.sams}") #| samtools view -Sb -F 4 -o {output.bams}
    #     else:
    #         shell("hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
    #         -1 {input.fq1} -2 {input.fq2} | samtools view -Sb -F 4 -o {output.bams}")
    shell:
        "hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
        -1 {input.fq1} -2 {input.fq2} | samtools view -Sb -F 4 -o {output.bams}"

# rule index_star:
#     input:
#         config["working_dir"] +"MorexV3/"+ config["ref"]
#     output:
#         [config["working_dir"]+ "MorexV3/" for i in [,,,,,]] #+ config["ref"].replace("fasta","") + str(i) + ".ht2" for i in range(1,9)]
#     message:
#         "indexing genome {output}"
#     params:
#         #index  = config["working_dir"] + "MorexV3/" + config["ref"].replace("fasta",""),
#         genomedir = config["working_dir"] + "MorexV3/",# + config["ref"]
#         genome = config["working_dir"] + "MorexV3/" + config["ref"],
#         gff = config["working_dir"] + "MorexV3/gene_annotation/" + config["gff3"]
#     threads: 1
#     shell:
#         "STAR --runThreadN 6 \
#         --runMode genomeGenerate \
#         --sjdbGTFtagExonParentTranscript Parent \
#         --genomeDir {params.genomedir} \
#         --genomeFastaFiles {params.genome} \
#         --sjdbGTFfile  {params.gff}\
#         --sjdbOverhang 99 \
#         --alignIntronMax 20000"
#
#
# rule run_STAR:
#     input:
#         fw = WORKING_DIR + "fastq/{sample}_1.fastq",
#         rev= WORKING_DIR + "fastq/{sample}_2.fastq",
#     output:
#         bam=protected(WORKING_DIR + "mapped/{sample}.bam"),
#         counts=WORKING_DIR + "mapped/{sample}.count.tab", #WORKING_DIR + "mapped/{sample}.bam"
#         log_file= "logs/star/{sample}.star.log",
#     params:
#         stranded= "--outSAMstrandField intronMotif",
#         prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
#         #outdir=lambda wildcards: "analysis/STAR/{sample}".format(sample=wildcards.sample),
#         readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
#         star_index = config["working_dir"]+ "MorexV3/",
#     threads: 2
#     message: "Running STAR Alignment on {wildcards.sample}"
#     benchmark:
#         "benchmarks/{sample}.run_STAR.txt"
#     shell:
#         "STAR --runMode alignReads --runThreadN {threads}"
#         " --genomeDir {params.star_index}"
#         " --readFilesIn {input} --readFilesCommand zcat"
#         " --outFileNamePrefix {params.prefix}."
#         " --outSAMstrandField intronMotif"
#         " --outSAMmode Full --outSAMattributes All {params.stranded}"
#         " --outSAMattrRGline {params.readgroup}"
#         " --outSAMtype BAM SortedByCoordinate"
#         " --limitBAMsortRAM 10000000000" #50 Gib
#         " --quantMode GeneCounts"
#         " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
#         " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"
# # rule sam2bam:
#     input:
#         sams  = WORKING_DIR + "mapped/{sample}.sam"
#     output:
#         bams  = WORKING_DIR + "mapped/{sample}.bam"
#     run:
#         shell("samtools view -Sb {input.sams} -F 4 -o {output.bams}")
# rule featureCounts:
# 	input:
# 		"star_align/{sample}/Aligned.sortedByCoord.out.bam"
# 	output:
# 		"featureCounts/{sample}"
# 	threads: 6
# 	log:
# 		"featureCounts/log/{sample}.log"
# 	params:
# 		R2_info
# 	shell:
# 		'''
# 		featureCounts.sh {input} {gdir} {threads} {params} {tmpdir}>&{log}
# 		'''
# rule featureCounts_all:
# 	input:
# 		expand("featureCounts/{sample}",sample=samples),
# 		"samples"
# 	output:
# 		"log/OK.featureCounts"
# 	log:
# 		"log/featureCounts.log"
# 	shell:
# 		'''
# 		cat samples | awk '{{print "featureCounts/"$0"\t"$0}}' >.samples_feature
# 		col_paste.sh  -s 1 .samples_feature >featureCounts.txt 2>{log}
# 		echo "Finished featureCounts" >{output}
