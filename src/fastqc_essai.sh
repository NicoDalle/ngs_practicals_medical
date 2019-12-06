#! /bin/bash

# Hsapiens
data="/home/rstudio/disk"
mkdir -p $data
cd $data
mkdir -p data_qc

cd data_qc

#Analyse qualité pour chaque fichier et enregistre les résultats dans data_qc

#Execute fastqc sur le fichier renseigné, enregistre le rapport dans data_qc
fastqc ../data_trimmed/paired/SRR3308950_trimmed_paired_1.fastq -o .
