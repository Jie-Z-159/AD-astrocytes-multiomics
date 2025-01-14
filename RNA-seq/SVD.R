#SVD for RNA-seq data

rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape2)
#preprocess for data ----
femData <-read.csv("gene_fpkm_A.csv",row.names = 1,header=TRUE,sep=",")
femData<-femData[,-c(31:39)]
#filter data
exp.mad<-apply(femData,1,mad)
selectgenes<-order(exp.mad,decreasing = T)[1:10000]
femData.filtererd<-femData[selectgenes,]
mydata<-femData.filtererd
#log
logmydata<-log10(mydata)
logmydata<-logmydata[is.finite(rowSums(logmydata)),]
#scale
z_scoredata<-scale(logmydata,center = TRUE, scale = TRUE)
z1<-z_scoredata[,2]
mean(z1)
sd(z1)
finitedata<-round(z_scoredata,5)
finitedata<-finitedata[is.finite(rowSums(finitedata)),]
#SVD----
svd_res<-svd(finitedata)
df_zscore.svd<-svd_res
df_zscore.svd.modes <- diag(df_zscore.svd$d)%*%t(df_zscore.svd$v)
df_singular_values <- data.frame(n=1:length(df_zscore.svd$d),singular_value=df_zscore.svd$d)

ggplot(df_singular_values,aes(x = n, y = singular_value)) + 
  geom_point(shape=1,color='#3299ff',size=1.5) +
  geom_line(color='#3299ff')+
  labs(y = 'Singular value', x = 'Singular value rank', title = 'Singular values of the astrocyte expression data')+
  theme_classic()+
  theme(plot.title=element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        text = element_text(size=12))

#mode1----
Modes_Comb <- df_zscore.svd.modes
nkeep <- 1
Time <- c(2, 4, 6, 9, 12)
# create empty matrix to hold t-test results
PT2Mat <- matrix(nrow=nkeep, ncol=5)
for (i in 1:nkeep) {
  
  XAD <- matrix(Modes_Comb[i,1:15], nrow=3)
  XWT <- matrix(Modes_Comb[i,16:ncol(Modes_Comb)], nrow=3)
  
  for (j in 1:5) {
    PT2Mat[i,j] <- t.test(XAD[,j], XWT[,j])$p.value
  }
  
  df <- data.frame(Time=Time, Group=rep(c('XAD', 'XWT'), each=5), Values=c(colMeans(XAD), colMeans(XWT)))
  
  p <- ggplot(df, aes(Time, Values, color=Group)) + 
    geom_point(size=2, shape=1) +
    geom_line(lwd=1) +
    labs(title=paste("Mode", i), x='month', y='Mean mode value') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major.y = element_line(linetype = "dashed", color = "gray"),
          panel.grid.minor = element_blank(),
          text = element_text(size=12)) +
    scale_color_manual(values=c("red", "black"), labels = c("AD", "WT")) +
    guides(colour=guide_legend(title="Group")) +
    scale_x_continuous(breaks = c(2, 4, 6, 9, 12))
  
max_value <- max(df$Values)
padding <- 0.003 * max_value  
y <- max_value - padding   

  # Add t-test results to the plot
  for (j in 1:5) {
    max_value <- max(df$Values)
    padding <- 0.001 * max_value  
    y <- max_value - padding    
    p_value <- PT2Mat[i, j]  
    if (p_value < 0.01) {
      label <- "**"  
    } else {
      label <- " "  
    }
    p <- p + annotate("text", x = Time[j], y = y, label = label, size = 5,color="red")
  }
  if (i == nkeep) {
    xlab("month")
  }
  print(p)
}

#first singular value corresponding sample loading----
V_matrix <- svd_res$v
first_singular_loading <- V_matrix[, 1]
group<-  sub("m.*", "m", colnames(femData.filtererd))
Sample <- factor(sub('A$', '',
                     colnames(femData.filtererd)),
                 levels = c("WT2m1", "WT2m2", "WT2m3", "AD2m1", "AD2m2", "AD2m3", 
                            "WT4m1", "WT4m2", "WT4m3", "AD4m1", "AD4m2", "AD4m3", 
                            "WT6m1", "WT6m2", "WT6m3", "AD6m1", "AD6m2", "AD6m3",
                            "WT9m1", "WT9m2", "WT9m3", "AD9m1", "AD9m2", "AD9m3",
                            "WT12m1", "WT12m2", "WT12m3", "AD12m1", "AD12m2", "AD12m3"))
df_first_singular_loading <- data.frame(Sample =Sample,
                                        Loading = first_singular_loading,
                                        group=group)
ggplot(df_first_singular_loading, aes(x = Sample, y = Loading,color=group)) +
  geom_point(size = 2) + 
  labs(y = ' Loading', x = 'Sample', 
       title = 'Mode1 loadings across samples') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        text = element_text(size=12),
        axis.text = element_text(angle = 45, hjust = 1,size = 10))+
  scale_color_manual(
    values = c("#8A2BE2", "#4682B4", "#FFD700", "#006400", 
               "#FF0000", "#DDA0DD", "#87CEEB", "#EEE8AA", "#98FB98", "#FFC0CB"))

#mode1 and total count
count<-read.csv(file="gene_count.csv",header = T,sep = ",",row.names = 1)
sequencing_depth <- colSums(count)
first_singular_loading <- V_matrix[match(colnames(count), colnames(femData.filtererd)), 1]
df_correlation <- data.frame(SequencingDepth = sequencing_depth, 
                             FirstSingularLoading = first_singular_loading)
cor<-cor(df_correlation$SequencingDepth,df_correlation$FirstSingularLoading,method = "pearson")

ggplot(df_correlation, aes(x = SequencingDepth, y = FirstSingularLoading))+
  annotate("text", x = max(df_correlation$SequencingDepth)*0.85 , y = max(df_correlation$FirstSingularLoading), 
           label = paste("Pearson correlation =", round(cor, 3)), size = 5, color = 'black')+
  geom_point(size = 2, color = '#3299ff') + 
  labs(y = 'Sample loading', x = 'Sequencing depth') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        text = element_text(size=12))

#scRNA top 10 genes in mouse model----
genes_of_interest <- c("PCBD1", "MIF", "ADIPOR2", "PRDX6", "EPHX1", "XIST", "G0S2", "HIBCH", "HSPH1", "HSP90AA1")
genes_of_interest<-data.frame(genes_of_interest)
topgenes<-read.csv(file = "/data/scRNAseq/top_genes.csv")
topgenes<-topgenes$Ensembl.Gene.ID
femData_genes <- transform(femData[rownames(femData) %in% topgenes, ], 
                           Gene = rownames(femData[rownames(femData) %in% topgenes, ]))
femData_long <- melt(femData_genes, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
#group
femData_long$Group <- ifelse(grepl("^AD", femData_long$Sample), "AD", "WT")
femData_long$Time <- as.numeric(sub(".*?(\\d+)m.*", "\\1", femData_long$Sample)) 
femData_summary <- aggregate(Expression ~ Gene + Time + Group, data = femData_long, FUN = mean)
p4 <- ggplot(femData_summary, aes(x = Time, y = Expression, color = Group, group = Group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  labs(y = 'Mean Expression Level', x = 'Time (months)',
       title = 'Expression of Selected Genes across Time Points') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)

#deseq2 result
df<-read.csv(file = "/data/RNAseq/AD vs WT/different month/6m/ADvsWT6m.csv",header = T,sep=",",row.names = 1)
df_data<-df[topgenes,]
df_data$Gene <- rownames(df_data)
df_data<-na.omit(df_data)
ggplot(df_data, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = padj < 0.05)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  coord_flip() + 
  labs(title = "log2FoldChange of Fatty Acid Metabolism Genes in AD vs WT",
       x = "Gene", y = "log2 Fold Change") +
  scale_fill_manual(name = "Significance", values = c("grey", "red"), labels = c("Not Significant", "Significant")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

