#!/bin/bash

#This script maps a number of files cotaining reads (fasta or fastq-format)
#To a Trinity assembly contig file. Then reports the abundance for each 
#contig and each read file. BWA is used for alignments

contigs=$1
reads=$2
threads=$3

bwa index $contigs

for read in $reads*.fastq; do
    echo $read
    stub=${read%.*}
    bwa mem -t $threads ${contigs} $read > ${stub}.sam
    samtools view -b -S ${stub}.sam > ${stub}.bam
    samtools sort ${stub}.bam -o ${stub}_sorted.bam 
    samtools index ${stub}_sorted.bam
    samtools idxstats ${stub}_sorted.bam > ${stub}_contig_abundances.txt 
done
