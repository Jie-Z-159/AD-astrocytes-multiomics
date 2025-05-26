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
### Pre-processing ----
mydata<-fqq
# Extract cell indices that meet the conditions
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

# Scatter plot
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
## Top 10 bar plot ----
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
    guide = guide_legend(override.aes = list(fill = NA, color = "#ff4856",stroke = 1))  # Modify legend circle to white
  ) +
  theme_minimal() +
  labs(x = "Metabolic pathway score", y = "Top 10 metabolic pathway") +  # Only capitalize first letter of y-axis label
  theme(
    text = element_text(size = 12, color = "black"),  # Global text size and color settings
    axis.title = element_text(size = 12, color = "black"),  # Axis titles
    axis.text = element_text(size = 12, color = "black"),  # Axis text
    legend.title = element_text(size = 12, color = "black"),  # Legend title
    legend.text = element_text(size = 12, color = "black"),  # Legend text
    legend.key = element_blank()
  )
print(p)

# Hallmark fatty acid ----
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
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +  # Add mean point
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.y = max(data_to_plot$MyModule1) * 1.1, vjust = -0.5) +  # Show significance symbol, adjust position
  annotate("text", x = 1.5, y = max(data_to_plot$MyModule1) * 1.2,
           label = "p-value = 8e-13", size = 5) +
  theme_minimal() +
  scale_fill_manual(values = c("#84b4e6", "#f2837f")) +  # Set fill colors, blue and red
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

# Change genes
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
    axis.text.x = element_text(angle = 45, hjust = 1), # Tilt X-axis text for better readability
    axis.title = element_text(size = 14),             # Adjust title font size
    plot.title = element_text(hjust = 0.5)            # Center plot title
  )

# Data description ----
DimPlot(Seurat, reduction = "umap",label = T, group.by = "Thal.phase",cols = c("Thal 0" = "#abdaf2", "Thal 2" = "#ea6163"))
meta_data <- Seurat@meta.data
donor_data <- meta_data %>%
  distinct(donor_id, .keep_all = TRUE) %>%  # Remove duplicate donor_id
  group_by(Thal.phase, sex) %>%                  # Group by ADNC and sex
  summarise(donor_count = n()) %>%         # Count donor numbers
  mutate(percentage = donor_count / sum(donor_count) * 100)  # Calculate percentage

sex_colors<-c("male"="#efcebb","female"="#a0cf8c")
ggplot(donor_data,aes(x=Thal.phase,y=donor_count,fill=sex ))+
  geom_bar(stat="identity", size=0.8)+
  coord_flip()+ # Make the plot horizontal
  theme_classic()+
  geom_text(aes(label=donor_count), vjust=0,hjust=-0.2, color="black", size=3.5)+
  theme(text = element_text(size = 14))+
  xlab("Number of donors")+
  ylab("Thal phase")+
  scale_fill_manual(values = sex_colors)

# Astrocyte marker ----
cell_marker<-read.csv(file="/data/scRNA allen/24-11/cell_marker.csv",header = T,sep = ",")
Seurat<-readRDS("Thal.RData")
cell_markers <- cell_marker %>% 
  group_by(Cell.name) %>% 
  summarise(markers = list(Cell.marker))

# Create an empty list to store expression results for each cell type
expression_results <- list()

# Loop through markers for each cell type
for (i in 1:nrow(cell_markers)) {
  cell_type <- cell_markers$Cell.name[i]
  markers <- cell_markers$markers[[i]]
  
  # Extract expression of these markers in single-cell data
  expression <- FetchData(Seurat, vars = markers,layer = "counts")
  
  # Calculate mean expression of these markers
  mean_expression <- rowMeans(expression, na.rm = TRUE)
  
  # Save results
  expression_results[[cell_type]] <- mean_expression
}

# Convert to data frame
expression_summary <- do.call(rbind, lapply(names(expression_results), function(cell_type) {
  data.frame(Cell.name = cell_type, MeanExpression = expression_results[[cell_type]])
}))

# Calculate average expression by cell type
summary_by_celltype <- expression_summary %>%
  group_by(Cell.name) %>%
  summarise(AverageExpression = mean(MeanExpression))

# Plot summary bar chart
summary_by_celltype$Cell.name <- factor(summary_by_celltype$Cell.name, 
                                        levels = summary_by_celltype$Cell.name[order(summary_by_celltype$AverageExpression)])

ggplot(summary_by_celltype, aes(x = Cell.name, y = AverageExpression, fill = Cell.name)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Expression of Cell Type Markers", x = "Cell Type", y = "Expression Level") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  coord_flip()
