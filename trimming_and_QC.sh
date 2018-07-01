#!/bin/sh

baseFileName=$1; # e.g. RJ-John-g-J_S2 /
startIndexR1=$2;
endIndexR1=$3;
startIndexR2=$4;
endIndexR2=$5;

# remove adapters and save only reads >= 23 bp (= 15bp + 8bp -> 4*2 random nucleotides), min quality=20 (both 5' and 3' ends)
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC --trim-n -q 20,20 -m 23 -o $baseFileName"_R1.trim1.fastq" -p $baseFileName"_R2.trim1.fastq" $baseFileName"_R1.fastq" $baseFileName"_R2.fastq" >> cutadapt.trim1.out 2>&1;

# trim poly-G tail
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -a "G{5}" -A "G{5}" -m 23 -O 5 -o $baseFileName"_R1.trim2.fastq" -p $baseFileName"_R2.trim2.fastq" $baseFileName"_R1.trim1.fastq" $baseFileName"_R2.trim1.fastq" >> cutadapt.trim2.out 2>&1;

# trim poly-A tail
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -a "A{5}" -A "A{5}" -m 23 -O 5 -o $baseFileName"_R1.trim3.fastq" -p $baseFileName"_R2.trim3.fastq" $baseFileName"_R1.trim2.fastq" $baseFileName"_R2.trim2.fastq" >> cutadapt.trim2.out 2>&1;

# trim adapters anew as adapters were not trimmed from poly-A-polyG tailed reads!
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -m 23 -O 5 -o $baseFileName"_R1.trim4.fastq" -p $baseFileName"_R2.trim4.fastq" $baseFileName"_R1.trim3.fastq" $baseFileName"_R2.trim3.fastq" >> cutadapt.trim4.out 2>&1;

#remove first and last 4 random nucleotides:
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -u 4 -u -4 -o $baseFileName"_R1.trim5.fastq" $baseFileName"_R1.trim4.fastq" >> cutadapt.trim5.1.out 2>&1;
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -u 4 -u -4 -o $baseFileName"_R2.trim5.fastq" $baseFileName"_R2.trim4.fastq" >> cutadapt.trim5.2.out 2>&1;

# verify length >= 15 paired data
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -m 15 -o $baseFileName"_R1.trim6.fastq" -p $baseFileName"_R2.trim6.fastq" $baseFileName"_R1.trim5.fastq" $baseFileName"_R2.trim5.fastq" >> cutadapt.trim6.out 2>&1;

# run fastqc and check if length-trimming is required
nohup /home/admin/software/FastQC/fastqc --outdir . $baseFileName"_R1.trim6.fastq" &> fastqc.trim6.R1.out;
nohup /home/admin/software/FastQC/fastqc --outdir . $baseFileName"_R2.trim6.fastq" &> fastqc.trim6.R2.out;

nohup /home/admin/software/fastx-toolkit/bin/fastx_trimmer -Q33 -f $startIndexR1 -l $endIndexR1 -i $baseFileName"_R1.trim6.fastq" -o $baseFileName"_R1.trim7.fastq";
nohup /home/admin/software/fastx-toolkit/bin/fastx_trimmer -Q33 -f $startIndexR2 -l $endIndexR2 -i $baseFileName"_R2.trim6.fastq" -o $baseFileName"_R2.trim7.fastq";

# re-run cutadapt to remove sequences that are too short (<15 bp)
nohup /home/admin/software/cutadapt-1.11/bin/cutadapt -m 15 -o $baseFileName"_R1.trimmed.fastq" -p $baseFileName"_R2.trimmed.fastq" $baseFileName"_R1.trim7.fastq" $baseFileName"_R2.trim7.fastq" >> cutadapt.trim8.out; 2>&1 &

# run final fastqc to check data
nohup /home/admin/software/FastQC/fastqc --outdir . $baseFileName"_R1.trimmed.fastq" &> fastqc.trimmed.R1.out;
nohup /home/admin/software/FastQC/fastqc --outdir . $baseFileName"_R2.trimmed.fastq" &> fastqc.trimmed.R2.out;

