#!/usr/bin/python3

import os
import sys
from subprocess import *

class quality_control():

    def __init__(self):
        self.fastq_dir = "/home/nazhang/luozhihui/project/CRC/RNA-seq/"
        self.Trimmomatic = "/home/nazhang/luozhihui/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
        self.outputDir = "/home/nazhang/luozhihui/project/CRC/trimmoResult"
        self.adaptor = "/home/nazhang/luozhihui/project/CRC/TruSeq2-PE.fa"

    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p

    def fastp(self, fastq1=None, fastq2=None, outputDir=None):
        sample = os.path.basename(fastq1).split(".")[0]
        output1 = sample + ".fastp.fq.gz"
        output1 = os.path.join(outputDir, output1)
        sample = os.path.basename(fastq2).split(".")[0]
        output2 = sample + ".fastp.fq.gz"
        output2 = os.path.join(outputDir, output2)
        sample = sample.split("_")[0]
        htmlPath = os.path.join(outputDir, sample + ".html")
        jsonPath = os.path.join(outputDir, sample + ".json")
        cmd = "fastp -z 4 -i %s -I %s -o %s -O %s -h %s -j %s" % (fastq1, fastq2, output1, output2, htmlPath, jsonPath)
        return cmd

    def extract(self, zipFile=None):
        """
        unzip the fastq file

        """

    def trimAdapterByTrimmomatic(self, fastq1=None, fastq2=None, outputDir=None):
        """
        fastq eg. "823-RNA_L8_1.fq.gz"
        ILLUMINACLIP:%s:2:30:10:
        %s is the adaptor file
        2 is bigest mismatch
        30 is palindrome method threshold value
        10表示simple方法的匹配阈值
        LEADING:3 表示切掉reads 5’端(the start of read)质量低于3的碱基或N
        TRAILING:3 表示切掉reads 3’端(the end of read)质量低于3的碱基或N
        SLIDINGWINDOW:4:15 表示以4个碱基作为窗口，窗口一个碱基一个碱基往后移，如果窗口内碱基的平均质量小于15,后面的都切掉
        MINLEN:36 #以上步骤处理后，如果reads的长度小于36，这条reads也会被排除

        :param fastq1: string
        :param fastq2: string
        :return:
        """
        fastq1_name = os.path.basename(fastq1)
        fastq2_name = os.path.basename(fastq2)
        forward_paired = os.path.join(outputDir, fastq1_name.replace("_1.fastp.fq.gz", "_1_paired.fq.gz"))
        forward_unpaired = os.path.join(outputDir, fastq1_name.replace("_1.fastp.fq.gz", "_1_unpaired.fq.gz"))
        reverse_paired = os.path.join(outputDir, fastq2_name.replace("_2.fastp.fq.gz", "_2_paired.fq.gz"))
        reverse_unpaired = os.path.join(outputDir, fastq2_name.replace("_2.fastp.fq.gz", "_2_unpaired.fq.gz"))
        logFile = os.path.join(outputDir, fastq1_name.replace("_1.fastp.fq.gz", ".log"))
        cmd = "java -jar %s PE -threads 20 -phred33 -trimlog %s\
         %s  %s %s %s %s %s \
         ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50" % \
              (self.Trimmomatic, logFile, fastq1, fastq2, forward_paired, forward_unpaired, reverse_paired, reverse_unpaired, self.adaptor)
        return cmd

    def fastqc(self, fastq=None, outputDir=None):
        #dirName = os.path.dirname(fastq)
        cmd = "fastqc -t 6 -o %s %s" % (outputDir, fastq )
        return cmd

    def star(self, fastq1=None, fastq2=None, outputDir=None):
        fastq1_name = os.path.basename(fastq1)
        prefix = os.path.join(outputDir, fastq1_name.replace("_1_paired.fq.gz", ""))
        cmd = "STAR --genomeDir star_genome --readFilesIn %s %s \
        --readFilesCommand zcat --outSAMstrandField intronMotif --runThreadN 8 --outFileNamePrefix %s" % \
              (fastq1, fastq2, prefix)
        return cmd

class count_read():

    def __init__(self):
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"


    def sam_to_bam(self, samFile=None, outputDir=None):
        fileName = os.path.basename(samFile)
        bamFile = os.path.join(outputDir, fileName.replace(".sam", ".bam"))
        cmd = "samtools view -@ 8 -F 4 -bS %s -o %s" % (samFile, bamFile)
        return cmd

    def sort_name(self, bamFile=None, outputDir=None):
        fileName = os.path.basename(bamFile)
        sortBamFile = os.path.join(outputDir, fileName.replace(".bam", ".sort.name.bam"))
        tmpDir = "/home/nazhang/luozhihui/project/CRC/_samtoolstmp"
        cmd = "samtools sort -n  -@ 8 -T %s -o %s %s" % (tmpDir, sortBamFile, bamFile)
        return cmd

    def htseq_count(self, sortBamFile=None, outputDir=None):
        fileName = os.path.basename(sortBamFile)
        resultTxt = os.path.join(outputDir, fileName.replace(".sort.name.bam", ".htseq.txt"))
        cmd = "htseq-count -f bam --stranded=no %s \
        /home/nazhang/luozhihui/project/CRC/reference_genome/Mus_musculus.GRCm38.94.gtf > %s" % (sortBamFile, resultTxt)
        return cmd

    def htseq_count_1(self, sortBamFile=None, outputDir=None):
        fileName = os.path.basename(sortBamFile)
        resultTxt = os.path.join(outputDir, fileName.replace(".sort.name.bam", ".htseq.txt"))
        cmd = "htseq-count -f bam --stranded=yes %s \
        /home/nazhang/luozhihui/project/CRC/reference_genome/Mus_musculus.GRCm38.94.gtf > %s" % (sortBamFile, resultTxt)
        return cmd

    def htseq_count_2(self, sortBamFile=None, outputDir=None):
        fileName = os.path.basename(sortBamFile)
        resultTxt = os.path.join(outputDir, fileName.replace(".sort.name.bam", ".htseq.txt"))
        cmd = "htseq-count -f bam --stranded=reverse %s \
        /home/nazhang/luozhihui/project/CRC/reference_genome/Mus_musculus.GRCm38.94.gtf > %s" % (sortBamFile, resultTxt)
        return cmd


if __name__ == "__main__":
    import re
    qc = quality_control()
    step = 8
    #step 1, generate fastqc file
    if step < 1:
        fastqFiles = os.listdir(qc.fastq_dir)
        for fi in fastqFiles:
            fastq = os.path.join(qc.fastq_dir, fi)
            cmd = qc.fastqc(fastq=fastq, outputDir="/home/nazhang/luozhihui/project/CRC/original_fastqc")
            qc.run(cmd=cmd)

    #step 2, run fastp
    if step < 2:
        fastqFiles = os.listdir(qc.fastq_dir)
        sample_num = []
        for file in fastqFiles:
            regx = re.compile("(\d+)-RNA_L8")
            result = regx.search(file)
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)
        for sample in sample_num:
            fastq1 = "%s-RNA_L8_1.fq.gz" % (sample)
            fastq1 = os.path.join(qc.fastq_dir, fastq1)
            fastq2 = "%s-RNA_L8_2.fq.gz" % (sample)
            fastq2 = os.path.join(qc.fastq_dir, fastq2)
            cmd = qc.fastp(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/fastp_step2")
            qc.run(cmd=cmd)

    #step 3, run trimomatic
    if step <3:
        fastq_dir = "/home/nazhang/luozhihui/project/CRC/fastp_step2"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)
        sample_num = []
        for file in fastqFiles:
            if re.search(".json", file):
                continue
            if re.search(".html", file):
                continue
            regx = re.compile("(\d+)-RNA_L8")
            result = regx.search(file)
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)
        #print(sample_num)
        #exit(1)
        for sample in sample_num:
            fastq1 = "%s-RNA_L8_1.fastp.fq.gz" % (sample)
            fastq1 = os.path.join(fastq_dir, fastq1)
            fastq2 = "%s-RNA_L8_2.fastp.fq.gz" % (sample)
            fastq2 = os.path.join(fastq_dir, fastq2)
            cmd = qc.trimAdapterByTrimmomatic(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/trimomatic_step3")
            qc.run(cmd)

    #step 4, run fastqc
    if step < 4:
        fastq_dir = "/home/nazhang/luozhihui/project/CRC/trimomatic_step3"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)
        for fi in fastqFiles:
            if re.search("fq.gz", fi):
                fastq = os.path.join(fastq_dir, fi)
                cmd = qc.fastqc(fastq=fastq, outputDir="/home/nazhang/luozhihui/project/CRC/fastqc_step4")
                qc.run(cmd=cmd)

    #step 5, run STAR
    if step < 5:
        fastq_dir = "/home/nazhang/luozhihui/project/CRC/trimomatic_step3"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)

        sample_num = []
        for file in fastqFiles:
            if re.search(".log", file):
                continue
            regx = re.compile("(\d+)-RNA_L8")
            result = regx.search(file)
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)

        for sample in sample_num:
            fastq1 = "%s-RNA_L8_1_paired.fq.gz" % (sample)
            fastq1 = os.path.join(fastq_dir, fastq1)
            fastq2 = "%s-RNA_L8_2_paired.fq.gz" % (sample)
            fastq2 = os.path.join(fastq_dir, fastq2)
            cmd = qc.star(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/star_step5")
            qc.run(cmd)

    bamPro = count_read()

    #step 6, convert sam to bam
    if step < 6:
        sam_dir = "/home/nazhang/luozhihui/project/CRC/star_step5"
        samFiles = os.listdir(sam_dir)
        if len(samFiles) == 0:
            exit(1)

        for file in samFiles:
            if not re.search(".sam", file):
                continue
            samfile = os.path.join(sam_dir, file)
            cmd = bamPro.sam_to_bam(samFile=samfile, outputDir="/home/nazhang/luozhihui/project/CRC/bam_step6")
            #print (cmd)
            qc.run(cmd)

    #step 7, sort bam
    if step < 7:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/bam_step6"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.sort_name(bamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/sort_step7")
            #print (cmd)
            qc.run(cmd)

    #step 8, htseq
    if step < 8:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step8")
            qc.run(cmd)

    if step < 9:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count_1(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step9")
            qc.run(cmd)

    if step < 10:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count_2(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step10")
            qc.run(cmd)
