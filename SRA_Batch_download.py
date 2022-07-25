import subprocess
#parallel -j 1 prefetch {} ::: $(cat SraAccList.txt)
import pandas as pd
import argparse
import os
import multiprocessing
import re

def pass_arg():
    parser = argparse.ArgumentParser(prog='SRA_Batch_download.py',description='Batch download SRA files for watermelon')
    parser.add_argument('--n', help='number of threads', type=int)
    parser.add_argument('--i', help='The SRA run file', type=str)
    parser.add_argument('--o', help='output folder for fastq', type=str)
    args = parser.parse_args()
    return args

# start multiple child and wget the data
def GET_sra(sra_file):
    output = "/home/yanjun/ct/data/sradata/"#args.o #+ sra_file.split("/")[-1]
    bashcmd = ["/home/yanjun/soft/sratoolkit.2.9.6-1-ubuntu64/bin/prefetch",sra_file,"-O", output]
    #bashcmd = ["wget",sra_file,output,"-O", output]
    cmd = subprocess.run(bashcmd,stdout=subprocess.PIPE,check=True)

def Batch_download(threads,sra_lists):
    p = multiprocessing.Pool(threads)
    p.map(GET_sra,sra_lists)

def sra2_fastq(inFile):
    sra_soft_path = "/home/yanjun/soft/sratoolkit.2.9.6-1-ubuntu64/bin/"
    output = "/home/yanjun/ct/data/sradata/fastq/"
    bashcmd = [sra_soft_path + "fastq-dump",inFile,"--split-files","-O", output,"--gzip"]
    cmd = subprocess.run(bashcmd,stdout=subprocess.PIPE,check=True)

def Batch_split(threads,inFiles):
    print( "number of threads is ", threads)
    p = multiprocessing.Pool(threads)
    p.map(sra2_fastq,inFiles)

def Get_SRA(File):
    #import pandas as pd
    df = pd.read_csv(File,sep="\t")
    df.sort_values(by="MBytes",inplace=True)
    Run = df["Run"].tolist()
    return Run[-200:-100]

def main():
    ''' Excute the work flow'''
    args = pass_arg()
    #print(args.i)
    sf = Get_SRA(args.i)
    #sfList = "https://sra-download.ncbi.nlm.nih.gov/traces/sra0/SRR/008546/" + sf
    Batch_download(threads=int(args.n),sra_lists =sf)
    print("File download completed")

def main_split():
    ''' Excute the work flow'''
    args = pass_arg()
    import re
    files = ["/home/yanjun/ct/data/sradata/"+f for f in os.listdir(args.o) if f.endswith(".sra")]
    #files = ["/home/yanjun/ct/data/sradata/"+f for f in os.listdir(args.o) if re.match('.*sra', f)]
    #in_Files =
    Batch_split(threads=args.n,inFiles=files)
    print("File file split completed")

#main
if __name__ == '__main__':
    #main()
    main_split()
hisat2 -p 2 --met-file /Users/yanjunzan/Documents/projects/barley_en/results/logs/ERR4393948_met.txt -x /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/Barley_MorexV3_pseudomolecules  -1 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/ERR4393948_R1_trimmed.fq.gz -2 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/ERR4393948_R2_trimmed.fq.gz | samtools view -Sb -F 4 -o /Users/yanjunzan/Documents/projects/barley_en/results/mapped/ERR4393948.bam
/Users/yanjunzan/OneDrive - Sveriges Lantbruksuniversitet/Projects/barley_en/bin/RNA_seq/Snakefile:
hisat2 -p 1 --met-file /Users/yanjunzan/Documents/projects/barley_en/results/logs/SRR8992621_met.txt -x /Users/yanjunzan/Documents/projects/barley_en/data/MorexV3/Barley_MorexV3_pseudomolecules             -1 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/SRR8992621_R1_trimmed.fq.gz -2 /Users/yanjunzan/Documents/projects/barley_en/results/trimmed/SRR8992621_R2_trimmed.fq.gz > /Users/yanjunzan/Documents/projects/barley_en/results/mapped/SRR8992621.sam
