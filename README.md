# MAGICAL

MAGICAL (Multiome Accessible Gene Integration Calling And Looping) analyzes paired scRNA-seq and scATAC-seq datasets from different conditions using a hierarchical Bayesian framework that improves model robustness by leveraging prior transcription factor motif and chromatin domain information. MAGICAL explicitly models errors both in chromatin accessibility and in gene expression measures, which improves the accuracy of the condition-specific 3D triads comprising distal chromatin sites, regulators, and their interactions with genes that are identified. 

Please check our paper "Mapping cell type regulatory triads modulated by disease from single-cell multiomics data" for more details. Any questions about the technical details of the paper or about this MAGICAL package can be emailed to: Xi Chen, xchen@flatironinstute.org

**MATLAB** and **R** scripts are provided. 


## Input files

Tools for scRNA-seq and scATAC-seq data processing are widely available, e.g. Seurat, ArchR. We realize that researchers may have different preference on data processing especially when there are multiple conditions, batches and samples involved. MAGICAL only needs gene symbols, peak coordinates, read count and cell meta information like cell type, sample/subject ID and sample group/condition. These information is very fundamental and should be easily obtained from any single cell multioimc dataset. We provide a ***R script*** to demo how to extra the following input files for MAGICAL from the intergated single cell data. 

**1. Cell type**

The scRNA-seq and scATAC-seq data sould be preprocessed and cell type labelled. MAGICAL infer regulatory triads for each cell type. Therefore, users need to specificy one cell type and then use the provided R script to prepare the following input files. 


**2. Candidate genes (DEG) and candidate chromatin sites (DAS)**:

As differential calling is usually done seperately during the scRNA-seq and scATAC-seq processing, we highly recommand users to prepare these two files using similar differential statistics cutoffs. For the selected cell type, the candidate gene file should be named as "(Cell type) candidate genes.txt" containing a list of gene symbols and the candidate chromatin site file should be named as "(Cell type) candidate regions.txt" including a list of regions (peaks) with genome coordinates (*chr, point1, point2*). We recommond running MAGICAL with hundreds of genes and a couple thousand of peaks. Too few genes like under 50 or two many genes like over 1000 will make the Bayesian process hard to converge.  

**2. scRNA-seq read count**

*Read count table*: a three-column matrix with *gene index*, *cell index*, and *RNA read count*

*Gene list*: a two-column matrix with 'gene index' and 'gene name'.

*Cell meta*: a five-column matrix with *cell index*, *cell barcode*, *cell type label*, *sample/subject ID*, and *condition* (can be more than two conditions but only up to two conditions will be analyzed in each run). Here, each sample must have a unique name and this name should be the same in the scATAC-seq data (to allow MAGICAL to pair data together). 


**3. scATAC-seq read count**

*Read count table*: a three-column matrix with *peak index*, *cell index*, and *ATAC read count*

*Peak list*: a four-column matrix with *peak index*, *chr*, *peak_point1*, *peak_point2*. The order of peaks (chr, point1, point2) should be corresponding to the order in the read count table

*Cell meta*: a five-column matrix with *cell index*, *cell barcode* (can be different from the scRNA-seq cell barcodes), *cell type label*, *sample/subject ID* (must be the same as scRNA-seq sample ID), and *condition* (must be the same as scRNA-seq condition)


**4. TF motif mapping (prior)**

*TF-peak mapping*: a three-column matrix with *peak index*, *motif index*, and *binary mapping*.

*Motif names*: a two-column matrix with *motif index* and *motif names*.


**5. TAD (prior)**

A six column matrix with genome coordinates (*left_chr, left_point1, left_point2, right_chr, right_point1, right_point2*) for the two boundaries of each domain. If no proper TAD information or HiC profile is available for the context being studied. We also provide another option to use distance to TSS as prior to initally pair peaks and genes. In addition, please ensure that the reference genome used for scATAC-seq and TAD is the same.  




## MAGICAL analysis

MAGICAL uses transcription factor (TF) motif and topological associated domains (TAD) as prior knowledge to infer regulatory triads of transcriptional regulators, regulatory chromatin sites and genes for each cell type. 

**Candidate triad constrcution:** To identify candidate disease-modulated triads, DAS within each cell type are associated with TFs by motif sequence matching. DAS are then linked to DEG in that cell type by requiring them to be within the same TAD or within a user controlled distance. During this step, to refine input genes and peaks, MAGICAL will run differential analysis using pseudobulk data. We set this step as optional but it has been domenstated to effectively lower the false positive rate in single cell differential calling (Nature Communication). 

**Triad linkage inference:** For each candidate triad, MAGICAL uses a Bayesian framework to iteratively model chromatin accessibility and gene expression variation across cells and samples in that cell type and estimate the strength of the triad TF-peak-gene linkages. The TF binding strength and TF activity are optimized to fit to the chromatin accessibility data. The estimated TF binding strength, TF activity and the gene expression data are used to infer the peak-gene interaction strength. We optimize the states of TF-peak-gene linkages based on the estimated strength which is used to initialize the next round of estimations. Finally, optimized triads fitting the variation in both data types are selected.

**Disease-associated triads output:** For each cell type, a file containing gene, chromatin site and regulator information will be finally produced by MAGICAL, with the name as "(Cell type) MAGICAL triads.txt". MAGICAL uses its default thresholds (posterior probabilities on TF-peak binding and peak-gene looping) to select triads and write them into the output file. Users can definitely adjust these thresholds in the provided scripts to allow more or fewer output triads. As the two linkages (TF-peak binding and peak-gene looping) in each triad are respectively identfied, we give higher priority on the peak-gene interaction when we select the final triads so it is possible to see some triads in the output file without high score TF bindings. These interactions are still important.  



