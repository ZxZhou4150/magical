# MAGICAL analysis

MAGICAL (Multiome Accessible Gene Integration Calling And Looping) is a hierarchical Bayesian approach that leverages paired scRNA-seq and scATAC-seq data from different conditions to map disease-associated transcription factors, chromatin sites, and genes as regulatory circuits. By simultaneously modeling signal variation across cells and conditions in both omics data types, MAGICAL achieved high accuracy on circuit inference. 

[R](https://github.com/Zxzhou4150/magical/tree/R-package) and [MATLAB](https://github.com/xichensf/magical/tree/main/MATLAB) scripts are provided for the use of MAGICAL with R v4.0.0 or MATLAB 2020a or later. The complete single cell data files are over 100GB and also the Bayesian learning process in MAGICAL may take hours to run on a local machine. 


![alt text](https://github.com/xichensf/magical/blob/main/MAGICAL.png)


## Demo
Before working on a complete dataset, we suggest users to download [our demo dataset](https://drive.google.com/file/d/1CerwMHMnS1PNFNMy00OoHQjn6T30M1j4/view?usp=sharing) and go through the demo tutorial for analysis ([R markdown](tutorial/MAGICAL.Rmd) and the [html page](https://zxzhou4150.github.io/MAGICAL.html)) and visualization ([R markdown](tutorial/Visualization.Rmd) and the [html page](https://zxzhou4150.github.io/Visualization.html)) to get familiar with MAGICAL functions (running time ~15 mins in total).


## Human PBMC single cell multiomics data
The sample-paired scRNA-seq and scATAC-seq data used in the MAGICAL paper have been deposited with the Gene Expression Omnibus under accession no. [GSE220190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220190). 

Users can also download datasets used in the papaer using the links below:
  * The Seurat-intergated R object of ***S. aureus*** scRNA-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/MRSA-MSSA-CTRL-all-combine-20210908.RData.gz) (15GB). 
  * The ArchR-integrated R project of ***S. aureus*** scATAC-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/Staph_scATAC_integration.tar.gz) (34GB, including arrow files).
  * The ArchR-integrated R project of **COVID-19** scATAC-seq data can be downloaded [here](https://wisp.princeton.edu/media/magical/COVID19_scATAC_integration.tar.gz) (7GB, including arrow files).


## Reference
Xi Chen et al., **Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data**, [*Nature Computational Science*, 3:644–657, (2023)](https://www.nature.com/articles/s43588-023-00476-5).



## Contact
Questions regarding the single cell multiome data and MAGICAL framework can be emailed to Xi Chen (<xchen@flatironinstitute.org>) and Yuan Wang (<yuanwang@cs.princeton.edu>). Other questions about our work should be emailed to Olga G. Troyanskaya (<ogt@genomics.princeton.edu>) and Stuart C. Sealfon (<stuart.sealfon@mssm.edu>).
