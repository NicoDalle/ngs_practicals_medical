#! /bin/bash

cd /home/rstudio/disk
mkdir -p star/index

#Liste des SRR d'intérêt
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

#Génération de l'index du génome humain annoté, avec 7 coeurs
#STAR --runThreadN 7 --runMode genomeGenerate \
#  --genomeDir star/index \
#  --genomeFastaFiles Hsap_genome.fa \
#  --sjdbGTFfile Hsap_annotation.gtf \
#  --sjdbOverhang 100

paired=/home/rstudio/disk/data_trimmed/paired

head $paired/SRR3308952_trimmed_paired_1.fastq

for srr in $SRR :
do
#Création d'un nouveau répertoire
mkdir star/$srr'_star'
cd star/$srr'_star'
#Quantification des reads
STAR --runThreadN 7 --genomeDir /home/rstudio/disk/star/index \
  --readFilesIn $paired/$srr'_trimmed_paired_1.fastq' \
  $paired/$srr'_trimmed_paired_2.fastq'
  
#Le fichier Aligned.out.sam est renvoyé par STAR, mais il est trop gros -> conversion en .bam, plus léger
samtools view -bS -h Aligned.out.sam > $srr'.bam'
#L'ancien fichier .sam est supprimé
rm Aligned.out.sam

done