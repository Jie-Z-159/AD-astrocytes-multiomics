#Use DEseq2 to find DEGs (example for 2m astrocytes)

rm(list=ls())
library(DESeq2)
#make count table
raw_df <- read.table(file = "gene_count.csv",header = T,sep = ",")
count_df <- raw_df[,c(5:7,20:22)]
head(count_df)
rownames(count_df) <- as.character(raw_df[,1])
colnames(count_df)
head(count_df)

# filter 
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df,1,function(x){ all(x > 0) }),]
dim(count_df)
dim(count_df.filter)
# condition table
sample_df <- read.table("condition2m.csv",header = T,sep = ",")
rownames(sample_df) <- c("AD2m1A" ,"AD2m2A", "AD2m3A", "WT2m1A", "WT2m2A", "WT2m3A")#改
print(sample_df)
sample_df$condition<-factor(sample_df$condition)
all(rownames(sample_df)==colnames(count_df.filter))

deseq2.obj <- DESeqDataSetFromMatrix(countData = count_df.filter, colData = sample_df, design = ~condition)
deseq2.obj
deseq2.obj$condition<-relevel(deseq2.obj$condition,ref = "WT")
deseq2.obj$condition
dds<-DESeq(deseq2.obj)
res<-results(dds)
res
res<-results(dds,contrast = c("condition","AD","WT"))
head(res)
resOrdered<-res[order(res$pvalue),]
resOrdered
sum(res$padj<0.1,na.rm = TRUE)
res0.05<-results(dds,alpha = 0.05)
summary(res0.05)
sum(res0.05$padj<0.05,na.rm = TRUE)
#MA plot
plotMA(res,ylim=c(-2,2))


vsd<-vst(dds,blind = FALSE)
rld<-rlog(dds,blind = FALSE)
ntd<-normTransform(dds)
head(assay(vsd),3)
head(assay(rld),3)
head(assay(ntd),3)

write.csv(as.data.frame(assay(vsd)),
          file = "count_transformation_VST.csv")
write.csv(as.data.frame(assay(rld)),
          file = "count_transformation_rlog.csv")

plotDispEsts(dds)
library(pheatmap)
select<-order(rowMeans(counts(dds,normalized=TRUE)),
              decreasing=TRUE)[1:100]
df<-as.data.frame(colData(dds)[,c("condition","type")])
mycolors=colorRampPalette(c("orange","white","purple"))(100)
pheatmap(assay(vsd)[select,],cluster_rows = FALSE,
         show_rownames = FALSE,
         border_color = NA,color = mycolors,cluster_cols = FALSE,
         annotation_col = df )
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition","type")])
mycolors = colorRampPalette(c("orange", "white","purple"))(100)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, 
         show_rownames=FALSE,
         border_color=NA,color=mycolors,cluster_cols=FALSE, annotation_col=df)
sampleDists<-dist(t(assay(vsd)))
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-paste(vsd$condition,vsd$type,sep = "-")
colnames(sampleDistMatrix)<-NULL
library("RColorBrewer")
colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,clustering_distance_cols = sampleDists,
         clustering_distance_rows = sampleDists,col=colors)
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
library(ggplot2)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

write.csv(as.data.frame(resOrdered),file = "ADvsWT2m.csv")
resSig<-subset(resOrdered,padj<0.05&abs(log2FoldChange) >= 0.5)
dim(resSig)
write.csv(as.data.frame(resSig),file = "dfADvsWT2m.csv")
resSigup<-subset(resOrdered,padj<0.05&log2FoldChange >= 0.5)
resSigdown<-subset(resOrdered,padj<0.05&log2FoldChange<= (-0.5))
dim(resSigup)
dim(resSigdown)
write.csv(as.data.frame(resSigup),file = "dfADvsWT2mup.csv")
write.csv(as.data.frame(resSigdown),file = "dfADvsWT2mdown.csv")
