# MAGICAL analysis

MAGICAL (Multiome Accessible Gene Integration Calling And Looping) analyzes paired scRNA-seq and scATAC-seq datasets from different conditions using a hierarchical Bayesian framework that improves model robustness by leveraging prior transcription factor motif and chromatin domain information. MAGICAL explicitly models errors both in chromatin accessibility and in gene expression measures, which improves the accuracy of the condition-specific 3D triads comprising distal chromatin sites, regulators, and their interactions with genes that are identified. 

Please check our paper "*Xi Chen et al. **Mapping disease-associated regulatory circuits by cell type from single-cell multiomics data**. 2022*" for more details. Any questions about the technical details of the paper or about this MAGICAL package can be emailed to: Xi Chen, xchen@flatironinstute.org.


## Input files

Tools for scRNA-seq and scATAC-seq data processing are widely available, e.g. Seurat, ArchR. We realize that researchers may have different preference on data processing especially when there are multiple conditions, batches and samples involved. MAGICAL only needs gene symbols, peak coordinates, read count and cell meta information like cell type, sample/subject ID and sample group/condition. These information is very fundamental and should be easily obtained from any single cell multioimc dataset. We provide a [R script](https://github.com/xichensf/magical/blob/main/Multiomics_input_for_MAGICAL.R) to demo how to extra the following input files for MAGICAL from the intergated single cell data. [(download demo input files)](https://drive.google.com/file/d/1CerwMHMnS1PNFNMy00OoHQjn6T30M1j4/view?usp=sharing)


#### **Cell type**

The scRNA-seq and scATAC-seq data sould be preprocessed and cell type labelled. MAGICAL infer regulatory triads for each cell type. Therefore, users need to specificy one cell type and then use the provided R script to prepare the following input files.


#### **Candidate genes (DEG) and candidate chromatin sites (DAS)**

As differential calling is usually done seperately during the scRNA-seq and scATAC-seq processing, we highly recommand preparing these two files using similar differential statistics cutoffs.  

  * *Candidate gene file*: a list of ``` gene symbols ```
  * *Candidate chromatin site file*: a three-column matrix of ```chr```, ```point1```, and ```point2``` 

#### **scRNA-seq read count**
We extract the scRNA-seq read count information from cells labelled to the selected cell type.   

  * *scRNA read count file*: a three-column matrix with ```gene index```, ```cell index```, and ```RNA read count```  
  * *scRNA gene name file*: a two-column matrix with ```gene index``` and ```gene name```.
  * *scRNA cell meta file*: a five-column matrix with ```cell index```, ```cell barcode```, ```cell type label```, ```sample/subject ID```, and ```condition```

Note, each sample must have a unique name and this name should be the same in the scATAC-seq data (to allow MAGICAL to pair data together). 


#### **scATAC-seq read count**
We extract the scATAC-seq read count information from cells labelled to the selected cell type. 

  * *scATAC read count file*: a three-column matrix with ```peak index```, ```cell index```, and ```ATAC read count```
  * *scATAC peak coordinate file*: a four-column matrix with ```peak index```, ```chr```, ```peak_point1```, ```peak_point2```.
  * *scATAC cell meta file*: a five-column matrix with ```cell index```, ```cell barcode```, ```cell type label```, ```sample/subject ID```, and ```condition```


#### **TF motif mapping (prior)**
We map 870 motifs from the [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs) library to all peaks and get their binary binding relationship. 

  * *Motif mapping file*: a three-column matrix with ```peak index```, ```motif index```, and ```binary binding```.
  * *Motif name file*: a two-column matrix with ```motif index``` and ```motif names```.

Note, the motif names of the same TF can be very different if they are collected from different databases. To avoid duplicated motifs, we require using the corresponding gene name for each TF motif. 

#### **TAD boundary (prior)**
Users will need to get the TAD boundary information from HiC profiles or similar experiments conducted in the same or similar context to their single cell datasets. We include a GM12878 cell line TAD file (~6000 domains with median size 400kb) in our demo for blood context analysis. 
  * *TAD file*: a three column matrix with ```chr```, ```left_boundary```, and ```right_boundary``` 

Alternatively, if no proper TAD information or HiC profile is available for the context being studied, we provide another option to use relative distance to TSS (e.g. 500kb) as prior to initally pair peaks and genes. Hg38 RefSeq file is provided for TSS reference.  


A demo:

```
loading all input data, it may take a while ...

We detected 2 conditions from the meta file.

The input scRNAseq data includes 36601 genes, 13248 cells from 12 samples/subecjts.

The input scATACseq data includes 144387 peaks, 13248 cells from 12 samples/subecjts.

There are paired data for 12 samples/subecjts. (check sample IDs if this number is lower than expected)

870 motifs, 11827 candidate chromatin sites and 1212 candidate genes are provided.
```


## Triad inference

MAGICAL uses TF motif and TAD as prior knowledge to infer regulatory triads of transcriptional regulators, regulatory chromatin sites and genes for the selected cell type. 

#### **1. Candidate triad constrcution**  
To identify candidate disease-modulated triads, candidate chromatin sites are associated with TFs by motif sequence matching. These sites are then linked to the candidate genes by requiring them to be within the same TAD or within a user controlled distance. 
```
Candidate regulatory circuits constrcution and MAGICAL model initialization

Candidate circuits include 81 TFs, 2379 chromatin sites, and 725 genes
```
#### **2. Triad linkage inference** 
For each candidate triad, MAGICAL uses a Bayesian framework to iteratively model chromatin accessibility and gene expression variation across cells and samples in that cell type and estimate the strength of the triad TF-peak-gene linkages. The TF binding strength and TF activity are optimized to fit to the chromatin accessibility data. The estimated TF binding strength, TF activity and the gene expression data are used to infer the peak-gene interaction strength. We optimize the states of TF-peak-gene linkages based on the estimated strength which is used to initialize the next round of estimations. Finally, optimized triads fitting the variation in both data types are selected.  
```
MAGICAL work starts ...

MAGICAL finished 10 percent

MAGICAL finished 20 percent

MAGICAL finished 30 percent

MAGICAL finished 40 percent

MAGICAL finished 50 percent

MAGICAL finished 60 percent

MAGICAL finished 70 percent

MAGICAL finished 80 percent

MAGICAL finished 90 percent

MAGICAL finished 100 percent
```
#### **3. Disease-associated triads output** 
For each cell type, a file containing gene, chromatin site and regulator information will be finally produced by MAGICAL, with the name as "(Cell type) MAGICAL triads.txt". MAGICAL uses its default thresholds (posterior probabilities on TF-peak binding and peak-gene looping) to select triads and write them into the output file. Users can adjust these thresholds in the provided scripts to allow more or fewer output triads. As the two linkages (TF-peak binding and peak-gene looping) in each triad are respectively identfied, we give higher priority on the peak-gene interaction when we select the final triads. Thus it is likely to see some triads in the output file without high score TF bindings. These interactions are still important.  



