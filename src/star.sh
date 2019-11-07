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

#génération de l'index du génome humain annoté, avec 7 coeurs
STAR --runThreadN 7 --runMode genomeGenerate \
  --genomeDir star/index \
  --genomeFastaFiles Hsap_genome.fa \
  --sjdbGTFfile Hsap_annotation.gtf \
  --sjdbOverhang 100

paired=/home/rstudio/disk/data_trimmed/paired

for srr in $SRR :
do
mkdir star/$srr'_star'
cd star/$srr'_star'
STAR --runThreadN 7 --genomeDir /home/rstudio/disk/star/index \
  --readFilesIn $paired/$srr'_trimmed_paired_1.fastq' \
  $paired/$srr'_trimmed_paired_2.fastq'
  
samtools view -bS -h Aligned.out.sam' > $srr'.bam'
rm Aligned.out.sam

done

