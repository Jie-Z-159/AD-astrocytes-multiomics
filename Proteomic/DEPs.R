#Use limma package to find DEPs

rm(list = ls())
library(ggplot2)
library(dplyr)
library(limma)
#load data
library(readxl)
rawdata<- read_xlsx("Protein quantitation.xlsx") 
ProteinIDs<-rawdata$`Protein IDs`
rawdata<-rawdata[,c(8:13)]
rownames(rawdata)<-ProteinIDs
write.csv(rawdata,file = "rawdata.csv")
data<-read.csv(file = "rawdata.csv",header = T,sep =",",row.names = 1)
#filter 
rm <- apply(data, 1, function(x){sum(x == 0) > 3})
df <- data[!rm,]
nrow(df)
#normalization
#log
df1 <- log2(df + 1)
boxplot(df1)
nor<-normalizeCyclicLoess(df1, weights = NULL, span=0.4, iterations = 3, method = "fast")
boxplot(nor)
write.csv(nor,file = "nor.csv")
#show intensity by boxplot
library(reshape2)
#ggplot
d1<-melt(df1)
p<-ggplot(d1,aes(x=as.factor(variable),y=value,fill=variable))
p+geom_boxplot()
d2<-melt(nor)
p<-ggplot(d2,aes(x=as.factor(Var2),y=value,fill=Var2))
p+geom_boxplot()


#group table
group<-read.csv(file = "group.csv",header = T,sep = ",",row.names = 1)
group_list <- factor(group[,1],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(nor)
cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
cont.matrix
#limma
fit <- lmFit(nor,design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
limma.na <- na.omit(tmpOut)
#DEPs
sum(limma.na$P.Value<0.05 & abs(limma.na$logFC) > 0.3)#
dif <- limma.na[limma.na$P.Value <= 0.05 & abs(limma.na$logFC) > 0.3,]#
down<- limma.na[limma.na$P.Value <= 0.05 & limma.na$logFC< -0.3,]#
up<-limma.na[limma.na$P.Value <= 0.05 & limma.na$logFC>0.3,]#

#output
write.csv(limma.na,file = "limma.na.csv")
write.csv(dif,file = "DEPs.csv")
write.csv(down,file = "down0.3.csv")
write.csv(up,file = "up0.3.csv")

#plot
library(ggplot2)
colnames(limma.na)
Dat<-limma.na
Dat$threshold = factor(ifelse(Dat$P.Value < 0.05 & abs(Dat$logFC) > 0.3, ifelse(Dat$logFC>0.3 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
ggplot(Dat,aes(x = logFC, y = -log10(P.Value),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(0.3,-0.3),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)

