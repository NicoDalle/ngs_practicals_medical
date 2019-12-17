# DESeq analysis Tp human
#DE Analysis:

# Libraries:
library("tximport")
library("readr")
library(apeglm)
library("DESeq2",quietly = T)
library("ggplot2")

# Location of the data:
dir <- "/home/rstudio/disk" ### Le dossier dans lequel vous travaillez

#Import des données 

#condition <- read.csv("~/disk/DESeq/condition.csv",header = T, sep=',') #on récupère les métadonnées
condition <- read.csv("~/disk/DESeq/metadata.csv",header = T, sep=',') #on récupère les métadonnées
samples <- data.frame(run=condition$Run,batch=condition$Batch,type=condition$subject_group,sex=condition$gender,time=condition$time,patient=condition$subject_id,GMB=condition$sample_from_primary_gbm_diagnosis)
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

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ type)

#Analyse des expressions différentielles

ddsTxi <- DESeq(ddsTxi)
res <-results(ddsTxi)
res

#Shrinkage des données
#On corrige l'effet du faible niveau d'expression sur le fold change

#On regarde le nom du coef qu'il faut mettre dans lfcShrink
resultsNames(ddsTxi)
#On fait le shrink
resLFC <- lfcShrink(ddsTxi, coef="type_responding_vs_non.responding", type="apeglm")
resLFC


#Classification des gènes selon leur p-value croissante

resOrdered <- res[order(resLFC$pvalue),]
resOrdered
summary(resLFC)

#Volcano plot
#On plotte le fold change en fonction du niveau d'expression
#Pour les gènes peu exprimés, les fold changes sont généralement grands, dû à des variations stochastiques
#Avec le shrinkage, on corrige pour cet effet

plotMA(res, ylim=c(-5,5))
plotMA(resLFC, ylim=c(-5,5))

#ACP

vsd <- vst(ddsTxi, blind=FALSE)
head(assay(vsd), 3)
#Plotter la PCA en basique, dans l'intgroup, on met tous les facteurs qu'on veut étudier
plotPCA(vsd, intgroup=c("sex", "type"))

#Plotter la PCA en joli
pcaData <- plotPCA(vsd, intgroup=c("batch", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

