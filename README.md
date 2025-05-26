### Multi-omics reveals changes in astrocyte fatty acid metabolism during early stages of Alzheimer's disease
We conducted an integrated multi-omics profiling of astrocytes obtained from APPswe/PSEN1ΔE9 transgenic AD and WT mice, 
including transcriptomics, proteomics, spatial metabolomics, to characterize the dynamic changes in astrocyte profiles over the course of AD progression.
To investigate whether similar changes are present in early human AD, we also analyzed single-nucleus RNA sequencing data of human brain samples.  

#### Bulk RNA-seq analysis
Differentially expressed genes (DEGs) between specified subgroups were identified utilizing the DESeq2 package.  
Singular value decomposition (SVD) analysis, gene expression levels were quantified using FPKM values and ranked based on the median absolute deviation (MAD) to identify the top 10,000 genes.  
Short Time-series Expression Miner (STEM)analysis, gene expression levels were quantified using FPKM values and ranked based on the MAD to identify the top 10,000 genes.  

#### Atrocyte-specific genome-scale metabolic model (GEM)
The astrocyte-specific GEM was constructed using the getINITModel function that implements the tINIT algorithm for reconstruction of cell-type-specific GEMs.

#### Proteomics analysis
The limma package was used for normalization and differential expression analysis between the WT and AD groups.

#### MALDI-MSI Data Analysis
We conducted MALDI-MSI on fresh frozen brain tissue from 6-month-old AD and WT mice with a resolution of 20 µm. 
The WT group (control) brain sections are WT1, AD2 (WT group), and WT3, while the AD group (APP/PS1 mice) brain sections are AD1, WT2 (AD group), and AD3.

#### snRNA-seq analysis
The snRNA-seq data utilized in this study were sourced from the Allen Institute for Brain Science and are available through 
the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) database (https://registry.opendata.aws/allen-sea-ad-atlas).  
For this study, data on astrocytes were obtained from 6 donors in the Thal 2 Aβ-deposition phase and 9 donors in the Thal 0 Aβ-deposition phase. All selected donors had no clinical diagnosis of dementia.

