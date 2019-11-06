#! /bin/bash

# Hsapiens
data="/home/rstudio/disk"
mkdir -p $data
cd $data
mkdir -p sra_data
cd sra_data

#Liste des SRR d'intérêt
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

for srr in $SRR
do
#Télécharge les séquences correspondant aux SRR renseignés avant
#split-file : séquençage pair-end, les séquences dans les 2 sens sont dans 2 fichiers différents
fastq-dump $srr -O /home/rstudio/disk/sra_data --split-files -I -X 100000

#Modifie les fichiers texte, NR: row number
#Change le nom des séquences en SRR..._1/1 et enlève la taille du read
awk -F"\." '{ if (NR%2 == 1) { $3="" ; print $1 "_" $2 "/1"} else {print $0} }' $srr'_1.fastq' > temp1.fastq
awk -F"\." '{ if (NR%2 == 0) { $3="" ; print $1 "_" $2 "/2"} else {print $0} }' $srr'_2.fastq' > temp2.fastq

#Enregistre les nouveaux fichiers à la place de anciens
mv temp1.fastq  $srr'_1.fastq'
mv temp2.fastq  $srr'_2.fastq'

done