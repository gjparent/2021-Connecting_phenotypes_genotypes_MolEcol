BiocManager::install("DESeq2")
BiocManager::install("DEFormats")
BiocManager::install("edgeR")

#Packages-----------------------
library("DESeq2")
library("DEFormats")
library("edgeR")
library("pheatmap")
library("ggplot2")
library("cluster")
library("factoextra")

#Data loading-------------------
rm(list=ls())
dev.off(dev.list()["RStudioGD"])

workingDir = "C:/Users/PARENT/Documents/R/2020-Adaptation-drought-white-spruce/Script/Results/";

data<- read.delim("C:/Users/PARENT/Documents/R/2020-Adaptation-drought-white-spruce/Data/Counts_PgDr.txt", row.names=1)
dataD <- read.delim("C:/Users/PARENT/Documents/R/2020-Adaptation-drought-white-spruce/Data/Meta_PgDr.txt", row.names=1)
countData <- as.matrix(read.delim("C:/Users/PARENT/Documents/R/2020-Adaptation-drought-white-spruce/Data/Counts_PgDr.txt", header=TRUE, row.names=1), nrow=37492,ncol=61)
colData <- read.delim("C:/Users/PARENT/Documents/R/2020-Adaptation-drought-white-spruce/Data/Meta_PgDr.txt", row.names = 1)
colData <- colData[,c("Clone", "Tx", "JJ")]
colnames(countData) <- NULL
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Clone + Tx + JJ, tidy=FALSE)
dds <- DESeqDataSet(dds, ~ Tx + JJ + Tx:Clone + Tx:JJ)
  
##Librairie sizes--------------
cl <- factor(dataD$Clone)
tx <- factor(dataD$Tx)
jj <- factor(dataD$JJ)
  
group <- factor(paste(dataD$Clone,dataD$Tx,dataD$JJ, sep="."))
design <- model.matrix(~0+group)
colnames(design)
d <- DGEList(counts=data,group=cl:tx:jj)
apply(d$counts, 2, sum) # total gene counts per sample
d$samples$lib.size <- colSums(d$counts)
LibSize <- d$samples
order(LibSize$lib.size, decreasing=TRUE)
barplot(d$samples$lib.size*1e-6, names=c(1:59), ylim=c(0,20), ylab="Library size (millions)")
#Export library size
write.csv(LibSize, file= "LibrariesSize_Drought.csv")

##Post-counts filtering--------------------------
dds <- dds[ rowSums(counts(dds)) > 1, ] #dim: 36701 59
dds <- dds[ rowSums(counts(dds)>5) >=2, ] #dim: 33824 59

##Estimate size factors and normalization--------
dds <- estimateSizeFactors(dds)
  
##LRT Model testing for the interaction between time and drought effect--------------------------------
ddsTC <- DESeqDataSet(dds, ~ Tx + JJ + Tx:JJ)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ Tx + JJ)
resTC <- results(ddsTC)
#Export DATA
resSig <- subset(resTC, padj < 0.05)
write.csv(as.data.frame(resSig), file="DESEQ_Dr_TxJJ.csv")
  
##Gene clustering of DEG Tx:JJ padj < 0.05 and 0.5 LFC-------------------------------------------------
  genelist1<-c(
    "PG_001495_T.1",
    "PG_006719_T.1",
    "PG_010167_T.1",
    "PG_000884_T.1",
    "PG_005555_T.1",
    "PG_001267_T.1",
    "PG_004784_T.1",
    "PG_011237_T.1",
    "PG_000110_T.1",
    "PG_014800_T.1",
    "PG_004766_T.1",
    "PG_001664_T.1",
    "PG_013894_T.1",
    "PG_024118_T.1",
    "PG_011865_T.1",
    "PG_010864_T.1",
    "PG_005136_T.1",
    "PG_013442_T.1",
    "PG_012858_T.1",
    "PG_011784_T.1",
    "PG_023211_T.1",
    "PG_005999_T.1",
    "PG_001171_T.1",
    "PG_011320_T.1",
    "PG_013253_T.1",
    "PG_005127_T.1",
    "PG_001676_T.1",
    "PG_013346_T.1",
    "PG_010373_T.1",
    "PG_001564_T.1",
    "PG_002005_T.1",
    "PG_008159_T.1",
    "PG_008276_T.1",
    "PG_014494_T.1",
    "PG_011200_T.1",
    "PG_013071_T.1",
    "PG_007303_T.1",
    "PG_001886_T.1",
    "PG_015270_T.1",
    "PG_001193_T.1",
    "PG_000584_T.1",
    "PG_015866_T.1",
    "PG_007685_T.1",
    "PG_005505_T.1",
    "PG_009724_T.1",
    "PG_001504_T.1",
    "PG_009729_T.1",
    "PG_006322_T.1",
    "PG_014168_T.1",
    "PG_006877_T.1",
    "PG_002362_T.1",
    "PG_010109_T.1",
    "PG_007295_T.1",
    "PG_002403_T.1",
    "PG_012672_T.1",
    "PG_004470_T.1",
    "PG_013905_T.1",
    "PG_010080_T.1",
    "PG_003313_T.1",
    "PG_014901_T.1",
    "PG_006047_T.1",
    "PG_002623_T.1",
    "PG_029236_T.1",
    "PG_014590_T.1",
    "PG_000630_T.1",
    "PG_014495_T.1",
    "PG_009337_T.1",
    "PG_013095_T.1",
    "PG_008056_T.1",
    "PG_022402_T.1",
    "PG_012815_T.1",
    "PG_001393_T.1",
    "PG_009553_T.1",
    "PG_012744_T.1",
    "PG_011998_T.1",
    "PG_007661_T.1",
    "PG_008028_T.1",
    "PG_014938_T.1",
    "PG_010927_T.1",
    "PG_014185_T.1",
    "PG_003105_T.1",
    "PG_009746_T.1",
    "PG_010849_T.1",
    "PG_015800_T.1",
    "PG_004484_T.1",
    "PG_010057_T.1",
    "PG_012475_T.1",
    "PG_005065_T.1",
    "PG_010143_T.1",
    "PG_011585_T.1",
    "PG_001877_T.1",
    "PG_014791_T.1",
    "PG_015723_T.1",
    "PG_006090_T.1",
    "PG_013630_T.1",
    "PG_016628_T.1",
    "PG_008434_T.1",
    "PG_000446_T.1",
    "PG_014761_T.1",
    "PG_012607_T.1",
    "PG_007453_T.1",
    "PG_014463_T.1",
    "PG_008461_T.1",
    "PG_006080_T.1",
    "PG_012439_T.1",
    "PG_010694_T.1",
    "PG_013057_T.1",
    "PG_005699_T.1",
    "PG_008285_T.1",
    "PG_006191_T.1"
  )

vsd <- varianceStabilizingTransformation(dds)
mat <- assay(vsd)[ genelist1, ]
mat <- mat - rowMeans(mat)
  
# Compute gap statistic
data1 <- mat
gap_stat <- clusGap(data1, FUN = pam,K.max = 50, B = 100)

# Plot gap statistic
fviz_gap_stat(gap_stat) # 2 clusters

#Heatmap
df <- as.data.frame(colData(vsd)[,c("Clone","Tx","JJ")])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
lut  = c(rev(cool), rev(warm))
annotation_colors = list(
    Clone = c("C08"="grey80","C11"="grey50", "C95"="black"),
    Tx = c(Control="dodgerblue2", Drought="yellow1"), JJ = c("D00"= "#9AFF9A","D07"= "#7BE27B","D14"= "#5DC45D","D18"= "#40A840","D22"= "#228B22"))
pl1<-pheatmap(mat, annotation_col=df,color = lut,border_color=FALSE, cluster_rows=TRUE, show_rownames=FALSE,show_colnames=FALSE,cluster_cols=FALSE,  clustering_distance_cols = "manhattan", clustering_method = "ward.D2", cellwidth = 6, cellheight = 2.5, annotation_colors = annotation_colors,fontsize_row = 4.5,cutree_rows=2, treeheight_row = 80)

hc <- pl1$tree_row
lbl <- cutree(hc, 2)
q <- data.frame(group = lbl, gene = names(lbl))
write.csv(as.data.frame(q),file="DESEQ_Dr_TxJJ_subset_groups.csv")