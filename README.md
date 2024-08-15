# multiomics-astrocyte-early-AD
This project includes R code for  transcriptomics, proteomics, spatial metabolomics, single-nucleus RNA sequencing data and NHANES data analysis. Scripts were implemented and tested under R version 4.1.1 (2021-08-10), R version 4.3.2 (2023-10-31). The R package vesion is,,,,
## RNA sequencing (RNA-seq) analysis
The RNA-seq dataset analyzed for this work was generated as part of another study, where it is described in detail and are accessible via the GEO series accession number GSE137028. The processed count table and FPKM table have been uploaded in the corresponding section. 
This part inclued three main analysis method: 
1. Differentially expressed genes (DEGs) between specified subgroups were identified utilizing the DESeq2. 
2. Tthe singular value decomposition (SVD) analysis.
3. the Short Time-series Expression Miner (STEM) analysis.
## Proteomics analysis
The raw proteomics analysis data could download from,,,.  The processed protein quantitation table has been uploaded in the corresponding section.
## matrix-assisted laser desorption/ionization mass spectrometry imaging (MALDI-MSI) data analysis
The raw MALDI-MSI data ccould download from,,,.  We had export the spot and spectra table for each brain slice. These data and identification feature table have been uploaded in the corresponding section.
## Single-nucleus RNA sequencing (snRNA-seq) analysis
The snRNA-seq data utilized in this study were sourced from the Allen Institute for Brain Science and are available through the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) database (https://registry.opendata.aws/allen-sea-ad-atlas)
## Analysis of dietary fatty acids and cognitive function
The dietary fatty acids and cognitive function data were download from the National Health and Nutrition Examination Survey (NHANES) 2011/2012, 2013/2014 two cycle.
