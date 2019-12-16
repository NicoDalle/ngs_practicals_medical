#! /bin/bash
echo "Need to add a dependency : "
#sudo apt-get install libxtst6 -y

cd /home/rstudio/disk
mkdir -p qualimap

#Liste des SRR d'intérêt (attention, on a ici que les SRR pas encore fait)
SRR="SRR3308952
SRR3308954
SRR3308956
SRR3308963
SRR3308971
SRR3308972
"

for srr in $SRR
do
cd /home/rstudio/disk/qualimap
mkdir -p $srr'_qualimap'

cd /home/rstudio/disk/star/$srr'_star'

#On trie le bam par ordre alphabétique des gènes, et pas selon leur position sur le génome
#-n : on trie par nom -O : type d'extension qu'on veut -o : output -@ : nombre de threads 
#dernier argument : fichier input
samtools sort -n -O bam -o $srr'_sorted.bam' -@ 7 $srr'.bam'

echo $srr

#Il faut trier les bam pour qualimap
#-bam : fichier bam d'input -gtf : fichier gtf d'annotation du génome -outdir : dossier de sortie
#-outfile : rapport de sortie -pe : paired end sequencing -s : le bam est trié (obligatoire comme -pe)
#qualimap rnaseq -bam $srr'_sorted.bam' \
#  -gtf /home/rstudio/disk/Hsap_annotation.gtf \
#  -outdir /home/rsudio/disk/qualimap/$srr'_qualimap' \
#  -outfile /home/rsudio/disk/qualimap/$srr'_qualimap'/$srr'_report.pdf' \
#  -pe -s

done

