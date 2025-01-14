#MSI data analysis

library(data.table)
library(Cardinal)
library(tibble)
library(dplyr)
library(tidyr)
library(Seurat)
library(sva)
library(mixOmics)
library(Seurat)
###proprecess-----
#example for AD1
#AD2, AD3,WT1,WT2, WT3 not shown
#load data 
AD1_spot<-read.csv(file = "AD1_spot.csv",sep = ";",skip=8)
AD_spot<-AD1_spot%>%select(-1)
AD1 <- fread("AD1-Spectra.csv",sep = ";",skip = 10,header = F)
AD<-AD1
#Transpose
AD1_t<- t(AD)
AD1_t[1:6,1:2]
AD1_t2<-AD1_t[-1,]%>%as.data.frame()%>%sapply(as.numeric)%>%as.data.frame()
#feature deduplication
AD1_t2<-arrange(AD1_t2,AD1_t2[[1]])
AD1_t3 <- AD1_t2 %>% distinct(V1, .keep_all = TRUE)
mz<-as.data.frame(AD1_t3[[1]])
which(duplicated(mz))
#pdata
coor<-AD_spot
coor <- data.frame(lapply(coor, as.integer))
run <- factor(rep("AD1", nrow(coor)))
diagnosis<-factor(rep("AD", nrow(coor)))
pdata <- PositionDataFrame(run=run, coord=coor,d=diagnosis)
# imagedata
intensit<-AD1_t3[,-1]%>%as.matrix()
colnames(intensit) <- NULL
str(intensit)
#featuredata
mz<- unlist(mz, use.names = FALSE)
fdata <- MassDataFrame(mz)
out1 <- MSImagingExperiment(imageData=intensit,
                            featureData=fdata,
                            pixelData=pdata)
pid <- pixels(out1,  y >= -13000, run == "AD1")
AD1_RAW<-out1[, pid]
image(AD1_RAW,mz=124.006)

###Combine different slices together----
load("rawdata.RData")
#Using feature identification table
feature_identification<-read.csv(file="feature identification.csv",header=T,sep=";")
merge<-arrange(feature_identification,m)
merge <- merge %>% distinct(m.z, .keep_all = TRUE)
mz<-as.data.frame(merge$m.z)
mz<- unlist(mz, use.names = FALSE)
raw<-cbind(AD1_RAW,AD2_RAW,AD3_RAW,WT1_RAW,WT2_RAW,WT3_RAW)
pro <- raw %>%peakBin(ref=mz,
                      tolerance=6.5,
                      units="ppm",
                      type="height") %>%
  process()

###import to seurat----
mse<-pro
intensity<-iData(mse)
# count
colnames(intensity) <- paste("V", 1:ncol(intensity), sep = "")
rownames(intensity)<-mse@featureData@mz
#metadata
run <- as.character(mse@elementMetadata@run)
d<-as.character(mse@elementMetadata@listData$d)
meta<-data.frame(run,d)
rownames(meta)<-paste("V",1:nrow(meta),sep="")
#object
se<-CreateSeuratObject(counts = intensity,project = "rms",meta.data = meta)
#add metadata
se@meta.data$plate <- ifelse(grepl("AD1|AD2|AD3", se@meta.data$run), "plate1",
                             ifelse(grepl("WT1|WT2|WT3", se@meta.data$run), "plate2", NA))
x<-mse@elementMetadata@coord$x
y<-mse@elementMetadata@coord$y
se@meta.data$x<-x
se@meta.data$y<-y
se$y <- ifelse(se$plate == "plate2", se$y - 10000, se$y)

###Swap Name----
log<-se
log@meta.data$temp <- log@meta.data$run
#Because all the original data have AD2 and WT2 reversed
# Replace AD2 with WT2
log@meta.data$run[log@meta.data$temp == "AD2"] <- "WT2"
# Replace WT2 with AD2
log@meta.data$run[log@meta.data$temp == "WT2"] <- "AD2"
# Delete temporary columns
log@meta.data$temp <- NULL

###normalization----
#Use 9AA as ref
reference_values <- FetchData(object = log, vars = "193.0758",layer="counts")
colnames(reference_values)<-"X193.0758"
reference_values <- reference_values %>%
  mutate(X193.0758 = ifelse(X193.0758 == 0, 1e-10, X193.0758))
reference_values <- as.numeric(reference_values$X193.0758)
# Get the expression matrix in the Seurat object
data_matrix <- log@assays[["RNA"]]@layers[["counts"]]
# Normalize each feature
normalized_matrix <- data_matrix / matrix(reference_values, nrow = nrow(data_matrix), ncol = ncol(data_matrix), byrow = TRUE)
#Manual tic
spotsum<-colSums(normalized_matrix)
ref_ticmatrix<-sweep(normalized_matrix,2,spotsum,FUN = "/")
ref_ticmatrix<-ref_ticmatrix*1000

###seurat preprocess----
# Store the normalized data back into the Seurat object
nornew<-CreateSeuratObject(ref_ticmatrix,meta.data = log@meta.data)
rownames(nornew)<-rownames(log)
summary(nornew@meta.data$nCount_RNA)
log<-NormalizeData(nornew,normalization.method = "LogNormalize",scale.factor =1000)
log <- FindVariableFeatures(log, selection.method = "vst", nfeatures = 444)
all.genes <- rownames(log)   
log <- ScaleData(log, features = all.genes)
log<- RunPCA(log,features =all.genes)
DimPlot(log, reduction = "pca",group.by="run",raster=FALSE)
explained_var <- log[["pca"]]@stdev^2 / sum(log[["pca"]]@stdev^2)
explained_var[1]
#umap setseed
log <- FindNeighbors(log, dims = 1:30)
log <- FindClusters(log, resolution = 0.1)
set.seed(2023082002)
log <- RunUMAP(log, dims = 1:30) 
DimPlot(log, reduction = "umap",raster = FALSE, group.by = "run", label = TRUE, label.size = 4)
colnames(log)<-colnames(se)
setwd("/data/MSI/MSI2/refticnew")
save(onlynor,file="onlynor.RData")

#combat remove batch effect----
data<-log@assays[["RNA"]]@layers[["data"]]
plate<-log@meta.data[["plate"]]
diagnose<-log@meta.data$d
design <- model.matrix(~ diagnose)
data_combat<-ComBat(dat=data,batch=plate,mod=design)
#data_combat<-ComBat(dat=data,batch=plate,mod=design,par.prior=F)
log_combat<-log
log_combat@assays[["RNA"]]@layers[["data"]] <- data_combat
log_combat <- FindVariableFeatures(log_combat, selection.method = "vst", nfeatures = 444)
all.genes <- rownames(log_combat)  
log_combat <- ScaleData(log_combat, features = all.genes)
log_combat<- RunPCA(log_combat,features =all.genes)
DimPlot(log_combat,reduction = "pca",group.by = "run",raster = F)
log_combat <- RunUMAP(log_combat, dims = 1:30) 
DimPlot(log_combat, reduction = "umap",raster = FALSE, group.by = "run", label = TRUE, label.size = 4)


###kmeans----
log<-log_combat
umap_data <- log[["umap"]]@cell.embeddings
centers<-5
kmeans_result <- kmeans(umap_data, centers=centers)
cluster_assignments <- kmeans_result$cluster
# Add clustering results to the Seurat object's metadata
log[["kmeansCluster"]] <- factor(cluster_assignments)
df <- data.frame(
  x =log@meta.data$x,
  y=log@meta.data$y,
  cluster = log@meta.data$kmeansCluster
)
#pdf("kmeans.pdf")
p <-ggplot(data = df,aes(x = x, y = y,fill=cluster)) +
  geom_tile(height=20,width=20)+theme_void() +
  ggtitle("kmeansCluster") +
  scale_fill_brewer(palette = "Spectral", type = "qual") +
  coord_fixed()
print(p)
#Combine 1 and 5 into 1
Idents(log) <- "kmeansCluster"
log<-RenameIdents(log,"5"="1")
levels(log)
log@meta.data[["kmeansCluster"]]<-log@active.ident
#dev.off()
combat_kmeans<-log
save(combat_kmeans,file="combat_kmeans.RData")

###ROI----
#example for AD1---- 
#AD3 and AD2（named by WT2）not shown
#Import the spot that circled before
file_list <- list.files(pattern = "_spots.csv$")
# Define a list to store data
data_list <- list()
# Loop through each file and process it
for (file in file_list) {
  data <- read.csv(file, skip = 8,sep=";") # Skip first 8 lines
  data_list[[file]] <- data
}
# Merge all data into one data frame
AD1roimerged_data<- do.call(rbind, data_list)
AD1roimerged_data <- data.frame(x = as.integer(AD1roimerged_data $x), y = as.integer(AD1roimerged_data$y))
###AD1 needs to be adjusted according to the coordinates of AD1
Idents(log)<-"run"
AD1<-subset(log,idents="AD1")
# Extract metadata from log object
AD1_metadata <- AD1@meta.data
AD1_metadata$colnames<-colnames(AD1)
# Use inner_join to combine log and AD1ROIspot based on x and y columns
merged_data <- inner_join(AD1_metadata, AD1roimerged_data, by = c("x", "y"))
# Use Seurat's subset function to obtain subsets
Idents(AD1)<-"kmeansCluster"
AD1_subset <- subset(AD1, idents="1")
#plot
df <- data.frame(
  x = AD1_subset@meta.data$x,
  y = AD1_subset@meta.data$y,
  ROI =AD1_subset@meta.data$region
)
p <-ggplot(data = df,aes(x = x, y = y,fill=ROI)) +
  geom_tile(height=20,width=20)+theme_void() +
  ggtitle("ROI") +
  scale_fill_brewer(palette = "Spectral", type = "qual") +
  coord_fixed()
print(p)
#Add region metadata to AD1
AD1$region<-"unselect"
AD1$region[AD1roimerged_data$colnames]<-"ROI"

##Combine the 3 slice ROIs together----
log_subset <- subset(log, cells =c(AD1roimerged_data$colnames,WT2roimerged_data$colnames,AD3roimerged_data$colnames))
#plot
df <- data.frame(
  x = log_subset@meta.data$x,
  y = log_subset@meta.data$y,
  nCount_RNA =log_subset@meta.data$nCount_RNA  
)
ggplot(data = df,aes(x = x, y = y,color=nCount_RNA)) +
  geom_tile(height=20,width=20)+theme_void() +
  labs("nCount_RNA")+ 
  scale_color_distiller(palette = "Spectral")+
  coord_fixed()
#Add region metadata to the log
log$region<-"unselect"
log$region[colnames(log_subset)]<-"ROI"
ROI<-log
save(ROI,file="ROI.RData")
###findmarkers for ROI
Idents(log)<-"kmeansCluster"
logs<-subset(log,idents="1")
Idents(logs)<-"region"
ROI_markers<-FindMarkers(logs,ident.1 = "ROI")
ROI_markers<-subset(ROI_markers,p_val<0.05)
ROI_markers<-merge(identi,ROI_markers,by="row.names",all=F)
ROI_markers<-drop_na(ROI_markers,Name)
ggplot(ROI_markers, aes(avg_log2FC, -log10(p_val))) + 
  geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, Name,"")), colour = "red",
                  size = 3,max.overlaps = Inf)

