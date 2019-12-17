# DESeq analysis Tp human
#DE Analysis:

# Libraries:
library("tximport")
library("readr")
library(apeglm)
library("DESeq2",quietly = T)

# Location of the data:
dir <- "/home/rstudio/disk" ### Le dossier dans lequel vous travaillez

#Import des données 

condition <- read.csv("~/disk/DESeq/condition.csv",header = T, sep=',') #on récupère les métadonnées
samples <- data.frame(run=condition$sample,type=condition$type,sex=condition$sex,time=condition$time,patient=condition$patient)
#file.path correspond au chemin pour arriver aux fichiers, comme avec des /
files <- file.path(dir,"salmon", samples$run, "quant.sf") #on va chercher tous les fichiers avec les données de salmon
names(files) <- samples$run

#On génère une matrice tx2gene qui fait le lien entre l'identifiant du transcrit et du gène correspondant

tx2gene <- as.character(read.table(files[1],header = T,sep = "\t")$Name)
trinity.genes <- unlist(lapply(lapply(strsplit(x = tx2gene,split = "|",fixed=T),FUN = `[`,1),paste,collapse="_"))
trinity.trans <- unlist(lapply(lapply(strsplit(x = tx2gene,split = "|",fixed=T),FUN = `[`,1:3),paste,collapse="|"))
tx2gene <- data.frame(txname=trinity.trans,geneid=trinity.genes)

#Génération de la table de compte
txi <- tximport(files,type="salmon",tx2gene=tx2gene)

#On peut afficher la table de compte avec :
txi$counts

ddsTxi <- DESeqDataSetFromTximport(txi, colData = condition, design = ~ type)
