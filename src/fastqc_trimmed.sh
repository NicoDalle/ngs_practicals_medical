#! /bin/bash

# Hsapiens
data="/home/rstudio/disk/data_trimmed/paired"
mkdir -p $data
cd $data
mkdir -p data_qc

#Récupére les noms des fichiers à analyser
FileList=`ls`

cd data_qc

#Analyse qualité pour chaque fichier et enregistre les résultats dans data_qc
for file in $FileList
do
echo $file
#Execute fastqc sur le fichier renseigné, enregistre le rapport dans data_qc
fastqc ../$file -o .
done