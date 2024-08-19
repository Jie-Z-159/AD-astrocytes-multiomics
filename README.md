# multiomics-astrocyte-early-AD
This project includes R code for  transcriptomics, proteomics, spatial metabolomics, single-nucleus RNA sequencing data and NHANES data analysis. Scripts were implemented and tested under R version 4.1.1 (2021-08-10), R version 4.3.2 (2023-10-31). The R package vesion is,,,,
## RNA sequencing (RNA-seq) analysis
The RNA-seq dataset was generated as part of another study, where it is described in detail and are accessible via the GEO series accession number GSE137028. The processed count table and FPKM table have been uploaded in the corresponding github section. 
This part inclued three main analysis method: 
1. Differentially expressed genes (DEGs) between specified subgroups were identified utilizing the DESeq2. 
2. Tthe singular value decomposition (SVD) analysis.
3. the Short Time-series Expression Miner (STEM) analysis.
## Proteomics analysis
The raw proteomics analysis data could download from iProX:IPX0009463000(https://www.iprox.cn/page/PSV023.html;?url=1724033062383bhPB).  The processed protein quantitation table has been uploaded in the corresponding github section.
## matrix-assisted laser desorption/ionization mass spectrometry imaging (MALDI-MSI) data analysis
The spot and spectra data could download from figshare(https://doi.org/10.6084/m9.figshare.26778319.v1). Specifically, six mouse brain sections were analyzed using MALDI-MSI. The WT group (control) brain sections are WT1, AD2 (WT group), and WT3, while the AD group (APP/PS1 mice) brain sections are AD1, WT2 (AD group), and AD3. 
Identification feature table and ROI table have been uploaded in the corresponding github section.
## Single-nucleus RNA sequencing (snRNA-seq) analysis
The snRNA-seq data utilized in this study were sourced from the Allen Institute for Brain Science and are available through the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) database (https://registry.opendata.aws/allen-sea-ad-atlas).
## Analysis of dietary fatty acids and cognitive function
The dietary fatty acids and cognitive function data were download from the National Health and Nutrition Examination Survey (NHANES) 2011/2012, 2013/2014 two cycle.
