#!/bin/sh

sampleName=$1; # e.g. RJ-John-g-J_S2

nohup blat Apis_mellifera.fa $sampleName"_R1.fasta" alignment/$sampleName"_R1.blast8.out" -q=dna -t=dna -out=blast8 &

nohup blat Apis_mellifera.fa $sampleName"_R2.fasta" alignment/$sampleName"_R2.blast8.out" -q=dna -t=dna -out=blast8 &


