#snRNA-seq analysis

library(Seurat)
packageVersion('Seurat')
library(dplyr)
library(patchwork)
#load data
snRNA<-readRDS("local.rds")
snRNA

#ensemble convert to gene name------
ref<-read.table(file="mart_export.txt",header = T,sep = ",")#download from ensembl website
library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
as_matrix()

class(snRNA@assays$RNA@counts)
dim(snRNA@assays$RNA@counts)
snRNA@assays$RNA@counts[1:5,1:5]
dense_matrix <- as_matrix(snRNA@assays$RNA@counts)
dense_matrix <- data.frame(dense_matrix)
dim(dense_matrix)
class(dense_matrix)
dense_matrix$ens <- row.names(dense_matrix)
head(dense_matrix$ens)
p <- match(rownames(dense_matrix),ref$Gene.stable.ID)
anno <- ref[p,]
dim(anno) 
dense_matrix$gene_symbol <- anno$Gene.name
dense_matrix$anno_ens <- anno$Gene.stable.ID
check = dense_matrix$anno_ens == row.names(dense_matrix)
which(check == FALSE)
which(check == TRUE)
dense_matrix_anno <- dense_matrix[!duplicated(dense_matrix$gene_symbol),]
dim(dense_matrix_anno)
which(is.na(dense_matrix_anno$"gene_symbol"))
which(is.numeric(dense_matrix_anno$"gene_symbol"))
dense_matrix_anno$"gene_symbol"[1500:1600]
which(dense_matrix_anno$"gene_symbol" == "") #Check the position where gene_symbol is an empty string
dense_matrix_anno[1588,70009:70012]
anno[which(anno$Gene.stable.ID == "ENSG00000083622"),] #View the situation of this part of the gene in the annotation file
which(is.na(dense_matrix_anno$gene_symbol))
dense_matrix_anno <- dense_matrix_anno[-which(is.na(dense_matrix_anno$gene_symbol)),] 
dense_matrix_anno <- dense_matrix_anno[-which(dense_matrix_anno$gene_symbol == ""),] 
dim(dense_matrix_anno)
rownames(dense_matrix_anno) <- dense_matrix_anno[,"gene_symbol"]
counts <- select(dense_matrix_anno,-c(ens,gene_symbol,anno_ens))
tail(colnames(dense_matrix_anno))
tail(colnames(counts))
dim(counts)
### Construct sample information metadata
metadt <- snRNA@meta.data
dim(metadt)
dim(counts)
colnames(counts) <- gsub('[.]', '-', colnames(counts))
#write.csv(counts,"counts.csv")
head(row.names(counts))
head(colnames(counts))
class(counts)
head(colnames(metadt))
head(row.names(metadt))
unique(metadt$assay_ontology_term_id)

### Constructing Seurat Objects-----
Seurat <- CreateSeuratObject(counts,meta.data = metadt)
#output counts
write.csv(counts,file="counts.csv")
##quality control
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#PercentageFeatureSet
Seurat [["percent.mt"]] <- PercentageFeatureSet(Seurat , pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(Seurat , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Seurat , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
Seurat  <- subset(Seurat , subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 4)
Seurat 
###normalization
Seurat <- NormalizeData(object = Seurat, normalization.method =  "LogNormalize",  scale.factor = 10000)
Seurat <- FindVariableFeatures(object = Seurat, 
                             selection.method = "vst", 
                             nfeatures = 2000)

scale.genes <- row.names(Seurat)
Seurat <- ScaleData(Seurat, features = scale.genes) 
### PCA
Seurat <- RunPCA(Seurat, features = VariableFeatures(Seurat)) 
pc.num = 1:30
### Finding neighboring data points
Seurat <- FindNeighbors(Seurat, dims = pc.num) 
### cluster
Seurat <- FindClusters(Seurat, resolution = 0.5) #
table(Seurat$seurat_clusters)
### UMAP
Seurat <- RunUMAP(Seurat, dims = pc.num)
DimPlot(Seurat, reduction = "umap",label = T, group.by = "Supertype")
### output
saveRDS(Seurat, file = "/data/scRNA allen/Seurat_fqq.rds")

#SCPA analysis-----
library(ComplexHeatmap)
library(SCPA)
library(Seurat)
library(tidyverse)
library(magrittr)
library(dyno)
library(ComplexHeatmap)
library(circlize)
#metabolic_pathways downlaod from SCPA package 
pathways2 <- "/data/scRNA allen/metabolic_pathways.csv"
tcm <- seurat_extract(Seurat,
                      meta1 = "ADNC", value_meta1 = "Low ADNC")
th1 <- seurat_extract(Seurat,
                      meta1 = "ADNC", value_meta1 = "Not AD")

set.seed(20240103)
metabolic_out <- compare_pathways(samples = list(tcm, th1), downsample = 5000,
                                   pathways = pathways2)

head(metabolic_out, 5)
metabolic_out<-read.csv(file = "metabolic_out.csv",header = T,sep = ",")
plot_rank(scpa_out = metabolic_out, 
          pathway = "HALLMARK_FATTY_ACID_METABOLISM", 
          base_point_size = 2, 
          highlight_point_size = 3)
pathway.reslut<-metabolic_out
#Scatter plot
best<-pathway.reslut[2,]
ggplot(pathway.reslut,aes(x = X, y =  qval)) + 
  geom_point(shape=1,color='#3299ff',size=1.5) +
  geom_line(color='#3299ff')+
  labs(y = 'Metabolic pathway score', x = 'Metabolic pathway rank', title = 'Single-cell metabolic pathway analysis')+
  theme_classic()+
  theme(plot.title=element_text(hjust = 0.5,size=15),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title = element_text(size=15))+
  geom_point(data=pathway.reslut[pathway.reslut$Pathway == "HALLMARK_FATTY_ACID_METABOLISM",],color="#f87b7e",size=3)+
  geom_label(  aes(label = "Hallmark fatty acid metabolism"), data=best,  nudge_y = 2,  alpha = 0.5,
               hjust=-0.05,vjust=3)

#Top10 Column Chart
pathway.reslut_10<-pathway.reslut[1:10,]
# Convert the text in the Pathway column to lowercase and replace underscores with spaces
pathway.reslut_10$Pathway <- gsub("_", " ", tolower(pathway.reslut_10$Pathway))
# Keep the values ​​in the adjPval column to three decimal places
options(digits = 3)
axis_text_colors <- rep("black", nrow(pathway.reslut_10))
axis_text_colors[9] <- "#df5975"  
axis_text_colors
ggplot(pathway.reslut_10, aes(y=qval, x=reorder(Pathway,(qval)) , fill=adjPval))+
  geom_bar(stat="identity")+
  scale_fill_continuous(low="#6ca0de", high='#b5ddf1',name="adj pval ")+
  theme_classic()+
  coord_flip() +
  labs(y="Metabolic pathway score",x="Top 10 metabolic pathway")+
  theme(text = element_text(size = 21),
        axis.text.y = element_text(color=axis_text_colors),
        axis.title.x = element_text(size = 19),
        axis.title.y=element_text(size = 19),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
# First, sort the data frame according to qval
pathway.reslut_10 <- pathway.reslut_10[order(pathway.reslut_10$qval, decreasing = F), ]
# Convert Pathway to a factor and specify the levels to maintain sort order
pathway.reslut_10$Pathway <- factor(pathway.reslut_10$Pathway, levels = pathway.reslut_10$Pathway)
# Making a ggplot bubble chart
p <- ggplot(pathway.reslut_10, aes(y = Pathway, x = qval, size = -log10(adjPval), fill = FC)) +
  geom_col() +
  geom_point(color = "#f19494",shape=21) +
  scale_fill_gradient(low = "#e9f3e5", high = "#a1d08d", name = "Fold Change (FC)") +
  scale_size_area(max_size = 13, name = "-log10(Adjusted P-value)") +
  theme_minimal() +
  labs(x = "Metabolic pathway score", y = "Top 10 metabolic pathway") +
  theme(
    text = element_text(size = 14, color = "black"), 
    axis.title = element_text(size = 16, color = "black"), 
    axis.text = element_text(size = 14, color = "black"), 
    legend.title = element_text(size = 16, color = "black"), 
    legend.text = element_text(size = 14, color = "black") 
  )

##Use a value to represent the fatty acid in each cell-----
fatty_acid_pathway<-read.csv(file="fatty acid metabolism.csv",header = T,sep = ",")
fattygenes<-fatty_acid_pathway$HALLMARK_FATTY_ACID_METABOLISM
subsub<-subset(Seurat,features=fattygenes,idents = c("Not AD", "Low"))
sum<-subsub@assays$RNA@scale.data%>%colSums()%>%as.data.frame()%>%rownames_to_column(var="cell")
submeta<-subsub@meta.data%>%as.data.frame()%>%select(,ADNC)%>%rownames_to_column(var="cell")
join<-inner_join(sum, submeta, by = "cell")
colnames(join)<-c("cell","scale_sum","ADNC")
join<-join%>%as.data.frame()%>%column_to_rownames(var="cell")
#boxplot
boxplot(scale_sum~ADNC,data=join,main="fatty acid metabolism")
p<-ggplot(data = join,mapping = aes(x=ADNC,y=scale_sum,fill = ADNC))+
  geom_boxplot()
#more information
get_box_stats <- function(y, upper_limit = max(join$scale_sum) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Cell =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n"
    )
  ))
}
p+stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  scale_fill_manual(values = c("#a5ddf8", "#faced2")) +
  theme_classic()+
  stat_summary(fun.data =get_box_stats, geom = "text", hjust = 0.5, vjust = 0.5)+
  ylab("Expression level(sacle data)")+
  theme(text = element_text(size = 14))+
  geom_signif(mapping=aes(x=ADNC,y=scale_sum), 
              comparisons = list(c("Not AD", "Low")), 
              map_signif_level=T,
              tip_length=c(0,0), 
              y_position = c(78), 
              size=1, 
              textsize = 4, 
              test = "t.test") 
  