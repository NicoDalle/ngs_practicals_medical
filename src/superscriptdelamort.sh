#! /bin/bash

# Hsapiens
data="/home/rstudio/disk"
mkdir -p $data
cd $data
mkdir -p sra_data_raw
mkdir -p sra_data
cd sra_data_raw

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
fastq-dump $srr -O /home/rstudio/disk/sra_data_raw --split-files -I -X 4

#Modifie les fichiers texte, NR: row number
#Change le nom des séquences en SRR..._1/1 et enlève la taille du read
awk -F"\." '{ if (NR%2 == 1) { $3="" ; print $1 "_" $2 "/1"} else {print $0} }' $srr'_1.fastq' > temp1.fastq
awk -F"\." '{ if (NR%2 == 1) { $3="" ; print $1 "_" $2 "/2"} else {print $0} }' $srr'_2.fastq' > temp2.fastq

#Enregistre les nouveaux fichiers à la place de anciens
mv temp1.fastq  ../sra_data/$srr'_1.fastq'
mv temp2.fastq  ../sra_data/$srr'_2.fastq'

done


cd /home/rstudio/disk

#Creation des dossiers pour stocker les fichiers nettoyés
mkdir -p data_trimmed
cd data_trimmed
mkdir -p paired
mkdir -p unpaired

cd ..

#Liste du numéro des patients
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

for fn in $SRR;
do
  #Trimmomatic sur des données Paired end
  	  #On utilise que 7 coeurs sur les 8 de la machine
  	  #Fichiers d'entrée
  	  #Fichiers de sortie
  	  #Fonctions à utiliser et arguments (les mêmes que dans le papier)
  	  
	java -jar /softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
		-threads 7 \
	  sra_data/$fn'_1.fastq' sra_data/$fn'_2.fastq' \
	  data_trimmed/paired/$fn'_trimmed_paired_1.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_1.fastq' \
	  data_trimmed/paired/$fn'_trimmed_paired_2.fastq' data_trimmed/unpaired/$fn'_trimmed_unpaired_2.fastq' \
	  ILLUMINACLIP:adaptateur.fa:2:30:10 \
	  LEADING:22 SLIDINGWINDOW:4:22 MINLEN:25

done

