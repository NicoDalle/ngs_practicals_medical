#! /bin/bash

cd /home/rstudio/disk/
mkdir -p star/index

#Liste des SRR d'intérêt
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

#Nettoyage des données du génome : on ne récupère que les chromosomes entiers
#Visualisation des scaffols non assemblés en chromosomes dans headers.txt
#grep ">" Hsap_genome.fa > headers.txt
#Récupération des headers des chromosomes entiers
#grep ">" Hsap_genome.fa |grep -v "_" |sed 's/>//g'>chr.txt
#Extraction des séquences des chromosomes assemblés en entier
#xargs samtools faidx Hsap_genome.fa < chr.txt> Hsap_chr.fa

#Le fichier Hsap_chr.fa est ensuite utilisé comme base pour l'index
#Il ne contient que les séquences qui sont bien assemblées



#Génération de l'index du génome humain annoté, avec 7 coeurs
STAR --runThreadN 7 --runMode genomeGenerate \
  --genomeDir star/index \
  --genomeFastaFiles Hsap_chr.fa \
  --sjdbGTFfile Hsap_annotation.gtf \
  --sjdbOverhang 100 \
  --genomeChrBinNbits 18


paired=/home/rstudio/disk/data_trimmed/paired

#for srr in $SRR :
#do
#Création d'un nouveau répertoire
#mkdir -p star/$srr'_star'
#cd star/$srr'_star'
#Quantification des reads
#STAR --runThreadN 7 --genomeDir /home/rstudio/disk/star/index \
#  --readFilesIn $paired/$srr'_trimmed_paired_1.fastq' \
#  $paired/$srr'_trimmed_paired_2.fastq'
  
#Le fichier Aligned.out.sam est renvoyé par STAR, mais il est trop gros -> conversion en .bam, plus léger
#samtools view -bS -h Aligned.out.sam > $srr'.bam'
#L'ancien fichier .sam est supprimé
#rm Aligned.out.sam

#done