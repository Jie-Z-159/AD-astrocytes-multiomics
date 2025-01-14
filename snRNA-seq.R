setwd("/data/scRNA allen/")
fqq<-readRDS("Seurat_fqq.rds")
library(Seurat)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(SCPA)
library(tidyverse)
library(magrittr)
library(circlize)
library(ggpubr)
###preprocess----
#load data
snRNA<-readRDS("local.rds") #download from SEA-AD
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
mydata<-Seurat

### Extract cells----
adnc_column <- mydata@meta.data$Cognitive.status
condition_indices <- which(adnc_column %in% c('No dementia'))
subset_data <- subset(mydata, cells = rownames(mydata@meta.data)[condition_indices])

adnc_column <- subset_data@meta.data$Thal.phase
condition_indices <- which(adnc_column %in% c('Thal 0', 'Thal 2'))
Seurat <- subset(subset_data,cells = rownames(subset_data@meta.data)[condition_indices])
saveRDS(Seurat,"Thal.RData")

DimPlot(Seurat, reduction = "umap",label = T, group.by = "Thal.phase")
Idents(Seurat)<-"Thal.phase"
table(Idents(Seurat))
th2markers<- FindMarkers(Seurat, ident.1 = "Thal 2",)
VlnPlot(Seurat, features = c("ADGRG7","FABP5"))
write.csv(th2markers,file="/data/scRNA allen/24-11/thal2markers.csv")
#GSEA----
library(org.Hs.eg.db)
deg<-th2markers
deg$symbol=rownames(deg)
df<-bitr(unique(deg$symbol),fromType="SYMBOL",
         toType=c("ENTREZID"),
         OrgDb=org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort<-DEG%>%
  arrange(desc(avg_log2FC))
geneList=data_all_sort$avg_log2FC
names(geneList)<-data_all_sort$ENTREZID
head(geneList)
str(geneList)
#kegg
kk2<-gseKEGG(geneList=geneList,
             organism='hsa',
             minGSSize=10,
             maxGSSize=200,
             pvalueCutoff=0.05)
kegg_result<-as.data.frame(kk2)
ridgeplot(kk2,10)
gseaplot2(kk2,
          title = "Lipid and atherosclerosis",  
          geneSetID = "hsa05417",               
          color = "red",                       
          base_size = 20,                       
          subplots = 1:2,                      
          pvalue_table = T)                    

#hallmark
names(geneList)<-data_all_sort$symbol
head(geneList)
str(geneList)

hm_H = msigdbr(species = "Homo sapiens",
               category = "H",
               subcategory = NULL) %>% as.data.frame() %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)
hm_H = hm_H %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
hm <- GSEA(geneList, TERM2GENE = hm_H)
hm<-as.data.frame(hm)



#SCPA----
pathways2 <- "/data/scRNA allen/metabolic_pathways.csv"
tcm <- seurat_extract(Seurat,
                      meta1 = "Thal.phase", value_meta1 = "Thal 2")
th1 <- seurat_extract(Seurat,
                      meta1 = "Thal.phase", value_meta1 = "Thal 0")

set.seed(20240103)
metabolic_out_new2024 <- compare_pathways(samples = list(tcm, th1), downsample = 4000,
                                          pathways = pathways2)

write.csv(metabolic_out_new2024,file="thal.csv")

plot_rank(scpa_out = pathway_result, 
          pathway = "HALLMARK_FATTY_ACID_METABOLISM", 
          base_point_size = 2, 
          highlight_point_size = 3)

astrocyte markers

####plot----
pathway_result<-read.csv(file = "thal.csv",header=T,sep=',')
pathway_result$FC<-pathway_result$FC * -1
colnames(pathway_result)[1] <- 'rank'
best <- pathway_result[pathway_result$Pathway == "HALLMARK_FATTY_ACID_METABOLISM",]
ggplot(pathway_result, aes(x = rank, y = qval)) + 
  geom_point(shape = 1, color = '#3299ff', size = 1.5) +
  geom_line(color = '#3299ff') +
  labs(
    y = 'Metabolic pathway score',
    x = 'Metabolic pathway rank',
    title = 'Single-cell metabolic pathway analysis'
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title = element_text(size = 15)
  ) +
  geom_point(
    data = best,
    color = "#f87b7e",
    size = 3
  ) +
  geom_label(
    data = best,
    aes(label = "Hallmark fatty acid metabolism"),
    nudge_y = 2,
    alpha = 0.5,
    hjust = -0.05,
    vjust = 3
  )
##top10 barplot----
pathway_result_10 <- pathway_result[order(pathway_result$qval, decreasing = TRUE), ][1:10, ]
pathway_result_10$Pathway <- str_to_title(gsub("_", " ", pathway_result_10$Pathway))
pathway_result_10$Pathway <- factor(pathway_result_10$Pathway, levels = rev(pathway_result_10$Pathway))
p <- ggplot(pathway_result_10, aes(y = Pathway, x = qval, size = -log10(adjPval), fill = FC)) +
  geom_col() +
  geom_point(color = "#ff4856", shape = 21) +
  scale_fill_gradient2(
    low = "#276419",
    mid = "#e5ffe5",
    high = "#8E0152",
    midpoint = 0,
    name = "Fold Change (FC)"
  ) +
  scale_size_area(
    max_size = 13,
    name = "-log10(Adjusted P-value)",
    guide = guide_legend(override.aes = list(fill = NA, color = "#ff4856",stroke = 1))  
  ) +
  theme_minimal() +
  labs(x = "Metabolic pathway score", y = "Top 10 metabolic pathway") +  
  theme(
    text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 12, color = "black"),  
    axis.text = element_text(size = 12, color = "black"),  
    legend.title = element_text(size = 12, color = "black"),  
    legend.text = element_text(size = 12, color = "black"),  
    legend.key = element_blank()
  )
print(p)

#hallmark fatty acid----
fatty_acid_pathway<-read.csv(file="fatty acid metabolism.csv",header = T,sep = ",")
gene_list<-fatty_acid_pathway$HALLMARK_FATTY_ACID_METABOLISM
Seurat<-AddModuleScore(
  object = Seurat,
  features = list(gene_list),
  name = "MyModule"
)
data_to_plot <- FetchData(Seurat, vars = c("MyModule1", "Thal.phase"))
p <- ggplot(data_to_plot, aes(x = Thal.phase, y = MyModule1, fill = Thal.phase)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.y = max(data_to_plot$MyModule1) * 1.1, vjust = -0.5) +  
  annotate("text", x = 1.5, y = max(data_to_plot$MyModule1) * 1.2,
           label = "p-value = 8e-13", size = 5) +
  theme_minimal() +
  scale_fill_manual(values = c("#84b4e6", "#f2837f")) +  
  labs(x = "Thal.phase", y = "Fatty acid pathway module Score") +
  theme(
    text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
print(p)
group_levels <- unique(data_to_plot$Thal.phase)
group1_data <- data_to_plot[data_to_plot$Thal.phase == group_levels[1], "MyModule1"]
group2_data <- data_to_plot[data_to_plot$Thal.phase == group_levels[2], "MyModule1"]

test_result <- wilcox.test(group1_data, group2_data)
print(test_result)

#change genes
Idents(Seurat) <- "Thal.phase"
gene_list_valid <- setdiff(gene_list, "CCDC58")
diff_exp_results <- FindMarkers(
  object = Seurat,
  ident.1 = "Thal 2",
  ident.2 = "Thal 0",
  features = gene_list_valid,
  logfc.threshold = 0,
  min.pct = 0
)%>%filter(p_val_adj<0.05)
#write.csv(diff_exp_results,file = "/data/scRNA allen/24-11/sctop_genes.csv")
RidgePlot(Seurat, features = c(  "HSP90AA1","HIBCH",
                               "MGLL", "GLUL"),ncol = 2) +
  scale_color_gradientn(colors = c("lightblue","lightpink")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),             
    plot.title = element_text(hjust = 0.5)            
  )


