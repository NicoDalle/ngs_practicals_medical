#! /bin/bash

# Hsapiens
data="/home/rstudio/disk"
mkdir -p $data
cd $data
mkdir -p data_qc

#Récupére les noms des fichiers à analyser
FileList=`ls sra_data/`

cd data_qc

#Analyse qualité pour chaque fichier et enregistre les résultats dans data_qc
for file in $FileList
do
echo $file
fastqc ../sra_data/$file -o .
done