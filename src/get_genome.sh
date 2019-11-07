cd /home/rstudio/disk

#Télécharge une URL, -O : output
wget -O Hsap_genome.fa.gz 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz'
wget -O Hsap_annotation.gtf.gz 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

#Décompression des fichiers .gz
gunzip Hsap_genome.fa.gz Hsap_genome.fa
gunzip Hsap_annotation.gtf.gz Hsap_annotation.gtf