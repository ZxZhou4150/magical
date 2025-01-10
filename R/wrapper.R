#' Main wrapper function
#' 
#' Main DAVINCI-MAGICAL wrapper function
#' 
#' @param RNA_counts RNA counts matrix. One row is one gene and one col is one spot.
#' @param ATAC_counts ATAC counts matrix. One row is one peak and one col is on spot.
#' @param niche_label The niche label array of DAVINCI output. Names must match colnames of the count matrices.
#' @param meta_add Additional meta data to add to the Seurat object. Rownames must match colnames of the count matrices.
#' @param pb Logic, pseudo-bulk level or not.
#' @param contrast Niche- or condition-specific contrast
#' @param niche1 (Required for all cases) For `contrast = "niche"`, this should be one niche you want to contrast. For `contrast = "condition"`, this should be the niche in which you want to contrast the conditions
#' @param niche2 The other niche you want to contrast if `contrast = "niche"`. Not specified is to contrast `niche1` with all other niches
#' @param condition The condition you want to make the contrast if `contrast = "condition"`
#' @param condition1 One condition you want to contrast if `contrast = "condition"`
#' @param condition2 The other condition you want to contrast if `contrast = "condition"`. Not specified is to contrast `condition1` with all other conditions.
#' @param p_thre Threshold of adjusted p-value (p_val_adj). Default is 0.05.
#' @param log2fc_thre Threshold of average log 2 fold change (avg_log2FC). Default is 0.3.
#' @param magical To perform downstream MAGICAL analysis or not. Default is `T`.
#' @param Ref_seq_file_path Path to the Refseq file for transcription starting site extraction
#' @param genome The genome for searching TF binding motifs. Default is `"hg38"`.
#' @param meta_spot_opt To run MAGICAL at meta-spot level or not. Default is `F`.
#' @param method How would you like to construct meta spots. Could be simple random, or based on the similarity of some features.
#' @param size The number of spots to be in one meta-spot
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' @param TAD_file_path The path to the TAD file.
#' @param dc Distance control parameter when TAD file is not specified, Default is 5e5 bps.
#' @param iteration_num The number of iterations for the main estimation step
#' @param Output_file_path The path to write the circuit results
#' @param ... Other parameters
#' 
#' @return A list of filtered data.frames in the form of the result of `Seurat::FindMarkers()`. The MAGICAL results are not returned but written into files.
wrapper_main = function(RNA_counts, ATAC_counts, niche_label, meta_add=NULL, pb, contrast = c("niche","condition"), niche1=NULL, niche2 = NULL, condition=NULL, condition1 = NULL, condition2 = NULL, p_thre = 0.05, log2fc_thre = 0.3, magical, Ref_seq_file_path, genome = "hg38", meta_spot_opt = F, method = c("simple_random", "feature"), size = 20, feature, TAD_file_path, dc = 5e5, iteration_num, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt', ...){
  ## step 1: differential analysis
  nsample = length(unique(obj$sample))
  if(nsample < 10){
    pb = F
    print("There are less than 10 samples. Spot level differential analysis will be applied.")
  }else{
    if(exists("pb") && (pb == F)){
      print("You are setting `pb = F`. There are more than 10 samples. Pseudo-bulk level differential analysis is recommended.")
    }else{
      pb = T  
      print("There are 10 or more samples. Pseudo-bulk level differential analysis will be applied. You can change to spot level by setting `pb = F`.")
    }
  }
  differentials = differential(RNA_counts, ATAC_counts, niche_label, meta_add, pb, contrast, niche1, niche2, condition, condition1, condition2, p_thre, log2fc_thre, ...)
  if(magical == T){
    ## step 2: prepare MAGICAL input
    loaded_data = prepare_magical_object(differentials[["deg"]], differentials[["das"]], RNA_counts, ATAC_counts, niche_label, meta_add, Ref_seq_file_path, meta_spot_opt, contrast, niche1, niche2, condition, condition1, condition2, method, size, feature, ...)
  
    ## step 3: run MAGICAL
    run_magical_main(TAD_file_path, dc, iteration_num, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt',...)
  }
  return(differentials)
}


#' Differential analysis
#' 
#' Differential gene and peak analysis with Seurat and Signac. Can be applied niche- or condition-specifically.
#' 
#' @param RNA_counts RNA counts matrix. One row is one gene and one col is one spot.
#' @param ATAC_counts ATAC counts matrix. One row is one peak and one col is on spot.
#' @param niche_label The niche label matrix of DAVINCI output. Rownames must match colnames of the count matrices.
#' @param meta_add Additional meta data to add to the Seurat object. Rownames must match that of `niche_label`.
#' @param pb Logic, pseudo-bulk level or not.
#' @param contrast Niche- or condition-specific contrast
#' @param niche1 (Required for all cases) For `contrast = "niche"`, this should be one niche you want to contrast. For `contrast = "condition"`, this should be the niche in which you want to contrast the conditions
#' @param niche2 The other niche you want to contrast if `contrast = "niche"`. Not specified is to contrast `niche1` with all other niches
#' @param condition The condition you want to make the contrast if `contrast = "condition"`
#' @param condition1 One condition you want to contrast if `contrast = "condition"`
#' @param condition2 The other condition you want to contrast if `contrast = "condition"`. Not specified is to contrast `condition1` with all other conditions.
#' @param p_thre Threshold of adjusted p-value (p_val_adj). Default is 0.05.
#' @param log2fc_thre Threshold of average log 2 fold change (avg_log2FC). Default is 0.3.
#' @param ... Other parameters for `Seurat::FindMarkers()`
#' 
#' @return A list of filtered data.frames in the form of the result of `Seurat::FindMarkers()` 
#'
# #' @import Seurat
# #' @import Signac
# #' @import dplyr
# #' @import tidyr
# #' @import purrr
#' 
#' @export
differential = function(RNA_counts, ATAC_counts, niche_label, meta_add=NULL, pb, contrast = c("niche","condition"), niche1=NULL, niche2 = NULL, condition=NULL, condition1 = NULL, condition2 = NULL, p_thre = 0.05, log2fc_thre = 0.3, ...){
  contrast = match.arg(contrast, c("niche","condition"))
  if(is.null(niche1))stop("Argument `niche1` is required for all cases.")
  
  #### check the input format
  #### test if metadata have those columns
  
  RNA_counts = RNA_counts[,names(niche_label)]
  ATAC_counts = ATAC_counts[,names(niche_label)]
  
  if(pb == T){
    ## build a Seurat object at pseudo-bulk level
    # aggregate counts
    spots_niche1id = which(niche_label == niche1)
    if(contrast =="niche"){
      RNA_niche1 = RNA_counts[,spots_niche1id]
      ATAC_niche1 = ATAC_counts[,spots_niche1id]
      meta_niche1 = meta_add[spots_niche1id,"sample"]
      if(is.null(niche2)){
        niche2 = paste0("non_",niche1)
        RNA_niche2 = RNA_counts[,-spots_niche1id]
        ATAC_niche2 = ATAC_counts[,-spots_niche1id]
        meta_niche2 = meta_add[-spots_niche1id,"sample"]
      }else{
        spots_niche2id = which(niche_label == niche2)
        RNA_niche2 = RNA_counts[,spots_niche2id]
        ATAC_niche2 = ATAC_counts[,spots_niche2id]
        meta_niche2 = meta_add[spots_niche2id,"sample"]
      }
      
      # aggregate to pb
      sample_df1 = data.frame(spot = names(meta_niche1), sample = meta_niche1)
      RNA_agg1 = aggregate_to_pb(sample_df1, RNA_niche1)
      ATAC_agg1 = aggregate_to_pb(sample_df1,ATAC_niche1)
      
      sample_df2 = data.frame(spot = names(meta_niche2), sample = meta_niche2)
      RNA_agg2 = aggregate_to_pb(sample_df2, RNA_niche2)
      ATAC_agg2 = aggregate_to_pb(sample_df2,ATAC_niche2)
      
      RNA = cbind(RNA_agg1, RNA_agg2)
      ATAC = cbind(ATAC_agg1, ATAC_agg2)
      
      # new metadata
      metadata = data.frame(
        spot = colnames(RNA),
        sample = c(unique(meta_niche1),unique(meta_niche2)),
        niche = rep(c(niche1,niche2),each = length(unique(meta_niche1)))
      )
      
      # build Seurat object
      obj <- Seurat::CreateSeuratObject(counts = RNA, meta.data = metadata)
      Seurat::DefaultAssay(obj) = "RNA"
      obj <- Seurat::NormalizeData(obj)
      
      atac_assay = Seurat::CreateAssayObject(counts = ATAC)
      obj[["ATAC"]] = atac_assay
      Seurat::DefaultAssay(obj) <- "ATAC"
      obj = Signac::RunTFIDF(obj)
      obj = Signac::FindTopFeatures(obj, min.cutoff = "q0")
      obj = Seurat::ScaleData(obj)
      
    }else if(contrast == "condition"){
      RNA_counts_totake = RNA_counts[,spots_niche1id]
      ATAC_counts_totake = ATAC_counts[,spots_niche1id]
      niche_label_totake = niche_label[spots_niche1id]
      meta_add_totake = meta_add[spots_niche1id,]
      
      spots_condition1id = which(meta_add_totake[[condition]] == condition1)
      
      RNA_condition1 = RNA_counts_totake[,spots_condition1id]
      ATAC_condition1 = ATAC_counts_totake[,spots_condition1id]
      meta_condition1 = meta_add[spots_condition1id,"sample"]
      if(is.null(condition2)){
        condition2 = paste0("non_",condition1)
        RNA_condition2 = RNA_counts_totake[,-spots_condition1id]
        ATAC_condition2 = ATAC_counts_totake[,-spots_condition1id]
        meta_condition2 = meta_add[-spots_condition1id,"sample"]
      }else{
        spots_condition2id = which(meta_add[[condition]] == condition2)
        RNA_condition2 = RNA_counts_totake[,spots_condition2id]
        ATAC_condition2 = ATAC_counts_totake[,spots_condition2id]
        meta_condition2 = meta_add[spots_condition2id,"sample"]
      }
      
      # aggregate to pb
      sample_df1 = data.frame(spot = names(meta_condition1), sample = meta_condition1)
      RNA_agg1 = aggregate_to_pb(sample_df1, RNA_condition1)
      ATAC_agg1 = aggregate_to_pb(sample_df1,ATAC_condition1)
      
      sample_df2 = data.frame(spot = names(meta_condition2), sample = meta_condition2)
      RNA_agg2 = aggregate_to_pb(sample_df2, RNA_condition2)
      ATAC_agg2 = aggregate_to_pb(sample_df2,ATAC_condition2)
      
      RNA = cbind(RNA_agg1, RNA_agg2)
      ATAC = cbind(ATAC_agg1, ATAC_agg2)
    }
  }else{
    # build a Seurat object at spot level
    metadata = cbind(niche_label, meta_add)
    
    # RNA
    obj <- Seurat::CreateSeuratObject(counts = RNA_counts, meta.data = metadata)
    Seurat::DefaultAssay(obj) = "RNA"
    obj <- Seurat::NormalizeData(obj)
    
    # ATAC
    atac_assay = Seurat::CreateAssayObject(counts = ATAC_counts)
    obj[["ATAC"]] = atac_assay
    Seurat::DefaultAssay(obj) <- "ATAC"
    obj = Signac::RunTFIDF(obj)
    obj = Signac::FindTopFeatures(obj, min.cutoff = "q0")
    obj = Seurat::ScaleData(obj)
  }
  
  if(contrast == "niche"){
    cells1 = obj$cell[which(obj$niche == niche1)]
    if(is.null(niche2)){
      cells2 = obj$cell[which(obj$niche != niche1)]
      }else{cells2 = obj$cell[which(obj$niche == niche2)]}
  }else if(contrast == "condition"){
    cells1 = obj$cell[which(obj$niche == niche1 & obj[[conditionname]] == condition1)]
    if(is.null(condition2)){
      cells2 = obj$cell[which(obj$niche == niche1 & obj[[conditionname]] != condition1)]
      }else{obj$cell[which(obj$niche == niche1 & obj[[conditionname]] == condition2)]}
  }
  
  diff_genes = Seurat::FindMarkers(obj@assays$RNA, cells.1 = cells1, cells.2 = cells2, ...)
  deg = diff_genes[which(diff_genes$p_val_adj<p_thre & abs(diff_genes$avg_log2FC)>logdc_thre),]
  
  diff_peaks = Seurat::FindMarkers(obj@assays$ATAC, cells.1 = cells1, cells.2 = cells2, ...)
  das = diff_peaks[which(diff_peakss$p_val_adj<p_thre & abs(diff_peaks$avg_log2FC)>logdc_thre),]
  
  return(list(genes = deg, peaks = das))
}

aggregate_to_pb = function(meta, mtx){
  agg_mtx = mtx %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(spot = rownames(.)) %>%
    left_join(meta, by = "spot") %>%
    select(-spot) %>%
    group_by(sample) %>%
    summarise(across(everything(), sum, .names = "{.col}")) %>%
    as.data.frame()
  rownames(agg_mtx) = agg_mtx$sample
  agg_mtx = t(as.matrix(agg_mtx[,-1]))
  return(agg_mtx)
}
 
#' Preparing MAGICAL input
#' 
#' Building MAGICAL object
#' 
#' @param deg Differentially expressed genes
#' @param das Differentially associated sites (peaks)
#' @param RNA_counts RNA counts matrix. One row is one gene and one col is one spot.
#' @param ATAC_counts ATAC counts matrix. One row is one peak and one col is on spot.
#' @param niche_label The niche label array of DAVINCI output. Names must match colnames of the count matrices.
#' @param meta_add Additional meta data to add to the Seurat object. Rownames must match colnames of the count matrices.
#' @param Ref_seq_file_path Path to the Refseq file for transcription starting site extraction
#' @param genome The genome for searching TF binding motifs. Default is `"hg38"`.
#' @param meta_spot_opt To run MAGICAL at meta-spot level or not. Default is `F`.
#' @param contrast Niche- or condition-specific contrast
#' @param niche1 (Required for all cases) For `contrast = "niche"`, this should be one niche you want to contrast. For `contrast = "condition"`, this should be the niche in which you want to contrast the conditions
#' @param niche2 The other niche you want to contrast if `contrast = "niche"`. Not specified is to contrast `niche1` with all other niches
#' @param conditionname The condition you want to make the contrast if `contrast = "condition"`
#' @param condition1 One condition you want to contrast if `contrast = "condition"`
#' @param condition2 The other condition you want to contrast if `contrast = "condition"`. Not specified is to contrast `condition1` with all other conditions.
#' @param method How would you like to construct meta spots. Could be simple random, or based on the similarity of some features.
#' @param size The number of spots to be in one meta-spot
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' @param ... Other parameters for `same_size_clustering()`
#' 
#' @import Matrix
#' @import Seurat
#' @import Signac
#' @import TFBSTools
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import chromVARmotifs
#'
#' @export
prepare_magical_object = function(deg, das, RNA_counts, ATAC_counts, niche_label, meta_add, Ref_seq_file_path, genome, meta_spot_opt = F, contrast = c("niche","condition"), niche1=NULL, niche2 = NULL, condition=NULL, condition1 = NULL, condition2 = NULL, method = c("simple_random", "feature"), size = 20, feature, ...){
  ## extract spots to use
  if(contrast == "niche" & !is.null(niche2)){
    totake = which(niche_label %in% c(niche1, niche2))
    RNA_counts = RNA_counts[,totake]
    ATAC_counts = ATAC_counts[,totake]
    niche_label = niche_label[totake]
    meta_add = meta_add[totake,]
  }
  else if(contrast == "condition"){
    if(is.null(condition2)){
      # take all the spots in `niche1`
      totake = which(niche_label == niche1)
      RNA_counts = RNA_counts[,totake]
      ATAC_counts = ATAC_counts[,totake]
      niche_label = niche_label[totake]
      meta_add = meta_add[totake,]
    }else{
      totake = which((niche_label %in% c(niche1, niche2)) & meta_add[[condition]] %in% c(condition1, condition2))
      RNA_counts = RNA_counts[,totake]
      ATAC_counts = ATAC_counts[,totake]
      niche_label = niche_label[totake]
      meta_add = meta_add[totake,]
    }
  }
  
  ## build MAGICAL object
  candidate_genes = data.frame(Gene_symbols = rownames(deg))
  
  candidiate_peaks_char = rownames(das)
  candidiate_peaks = do.call(rbind, strsplit(candidiate_peaks_char, "-"))
  candidiate_peaks = data.frame(chr = candidiate_peaks[,1], point1 = as.numeric(candidiate_peaks[,2], point2 = as.numeric(candidate_peaks[,3])))
  
  RNA_count_mtx = as(RNA_counts, "TsparseMatrix")
  
  RNA_genes = data.frame(Gene_index = 1:nrow(RNA_counts), Gene_symbols = rownames(RNA_counts))
  
  ATAC_count_mtx = as(ATAC_counts, "TsparseMatrix")
  
  ATAC_peaks_char = rownames(ATAC_counts)
  ATAC_peaks = do.call(rbind, strsplit(ATAC_peaks_char, "-"))
  ATAC_peaks = data.frame(chr = ATAC_peaks[,1], point1 = as.numeric(ATAC_peaks[,2], point2 = as.numeric(candidate_peaks[,3])))
  
  Ref_seq = read.table(Ref_seq_file_path, header = TRUE, sep = "\t")
  colnames(Ref_seq) = c("chr", "strand", "start", "end", "Gene_symbols")
  
  # get motifs with chromVARmotifs
  chromatinassay <- CreateChromatinAssay(counts = ATAC_counts, genome = genome)
  object <- CreateSeuratObject(counts = chromatinassay)
  object <- AddMotifs(object = object, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = human_pwms_v2)
  Peak_motif_mapping <- as(object@assays$ATAC@motifs@data * 1, "TsparseMatrix")
  motif_mapping <- summary(Peak_motif_mapping)
  motifs <- Peak_motif_mapping@Dimnames[[2]]
  motifs <- sapply(motifs, function(x) {
    splt <- strsplit(x, split = "_")
    if (grepl("ENSG", splt[[1]][1])) {
      splt[[1]][3]
    } else {
      splt[[1]][4]
    }
  }) # convert to gene names
  motifs <- cbind.data.frame(seq_len(length(motifs)), motifs)
  colnames(Motifs) <- c("motif_index", "name")
  
  ## Set metadata (RNA_cell and ATAC_cells) as indicated by meta_spot_opt
  if(meta_spot_opt == F){
    if(contrast == "niche"){
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche_label, subject_ID = meta_add$sample, condition = "1")
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche_label, subject_ID = meta_add$sample, condition = "1")
    }else if(contrast == "condition"){
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche1, subject_ID = meta_add$sample, condition = meta_add[[condition]])
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche1, subject_ID = meta_add$sample, condition = meta_add[[condition]])
    }
  }else{
    # labels (niche or condition)
    if(contrast == "niche"){
      clusters = niche_label
      if(is.null(niche2)){cluster[which(cluster!=niche1)] = paste0("non_",niche1)}
    }else if(contrast == "condition"){
      clusters = meta_add[[condition]]
      if(is.null(condition2)){cluster[which(cluster!=condition1)] = paste0("non_",condition1)}
    }
    
    method = match.arg(method, c("simple_random", "feature"))
    
    if(method == "simple_random"){
      new_idents = metaspot_sr(clusters, size)
    }else{
      new_idents = metaspot_feature(clusters, size, feature)
    }
    if(contrast == "niche"){
      #### hierarchy: samplle, niche??
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche_label, subject_ID = new_ident, condition = "1")
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche_label, subject_ID = new_ident, condition = "1")
    }else if(contrast == "condition"){
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche1, subject_ID = new_ident, condition = meta_add[[condition]])
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche1, subject_ID = new_ident, condition = meta_add[[condition]])
    }
  }
  
  Common_samples <- intersect(RNA_cells$subject_ID, ATAC_cells$subject_ID)
  
  loaded_data = list(
    "Common_samples" = Common_samples,
    "Candidate_Genes" = data.frame(candidate_Genes),
    "Candidate_Peaks" = data.frame(candidate_Peaks),
    "scRNA_Genes" = data.frame(RNA_genes),
    "scRNA_cells" = data.frame(RNA_cells),
    "scRNA_read_count_matrix" = RNA_count_mtx,
    "scATAC_Peaks" = data.frame(ATAC_peaks),
    "scATAC_cells" = data.frame(ATAC_cells),
    "scATAC_read_count_matrix" = ATAC_count_mtx,
    "Motifs" = data.frame(Motifs),
    "TF_Peak_binding_matrix" = motif_mapping,
    "Refseq" = data.frame(Refseq)
  )
  
  return(loaded_data)
}

#' Functions for metaspot construction
#'
#' Simple random
#' 
#' @seealso [metaspot_feature()]
#'
#' @param clusters The labels of the spots
#' @param size The size you want each metaspot group to have
#' 
#' @return A `data.frame` with 3 columns: index, label (niche or condition), group_index (assigned metaspot group)
#' 
#' @export
metaspot_sr <- function(clusters, size) {
  results <- data.frame(index = integer(), label = character(), group_index = integer())
  labels <- unique(clusters)
  
  for (label in labels) {
    group_number <- 1
    totake <- which(clusters == label)
    if(length(totake)==0){
      next()
    }
    
    subset_clusters <- cbind.data.frame(1:length(totake), clusters[totake])
    
    # randomize the clusters
    subset_clusters <- subset_clusters[sample(nrow(subset_clusters)), , drop = FALSE]
    
    num_items <- nrow(subset_clusters)
    num_full_groups <- num_items %/% size
    remainder <- num_items %% size
    
    for (i in seq_len(num_full_groups)) {
      indices <- (i - 1) * size + 1:size
      subset_clusters[indices, "group_index"] <- rep(paste0(label,"-",group_number), length(indices))
      group_number <- group_number + 1
    }
    
    if (remainder > 0) {
      indices <- (num_full_groups * size + 1):num_items
      if (remainder <= size / 3 & group_number != 1) {
        subset_clusters[indices, "group_index"] <- rep(paste0(label,"-",group_number - 1), length(indices))
      } else {
        subset_clusters[indices, "group_index"] <- rep(paste0(label,"-",group_number), length(indices))
      }
    }
    
    results <- rbind(results, subset_clusters)
  }
  return(results)
}

#' Functions for metaspot construction
#'
#' Based on the similarity of features (e.g., LVs, coordinates)
#' 
#' @seealso [metaspot_sr()]
#' 
#' @param clusters The labels of the spots
#' @param size The size you want each metaspot group to have
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' 
#' @return A `data.frame` with 3 columns: index, label (niche or condition), group_index (assigned metaspot group)
#' 
#' @export
metaspot_feature <- function(clusters, size, feature) {
  results <- data.frame(index = integer(), label = character(), group_index = integer())
  labels <- unique(clusters)
  
  for (label in labels) {
    totake <- which(clusters == label)
    if(length(totake)==0){
      next()
    }
    
    subset_clusters <- cbind.data.frame(1:length(totake), clusters[totake])
    
    if(length(totake)<=size){
      subset_clusters$group_index <- "1"
    }
    else{
      subset_features <- feature[subset_clusters[,1], ]
      result <- same_size_clustering(subset_features, clsize = 10)
      subset_clusters$group_index <- paste0(label, "-", result)
    }
    results <- rbind(results, subset_clusters)
  }
  return(results)
}

## functions from `same_size_clustering.R`
# https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd

#' Same Size Clustering
#'
#' This is a wrapper for several implementation that classify samples into
#' same size clusters, the details please see [this blog](http://jmonlong.github.io/Hippocamplus/2018/06/09/cluster-same-size/).
#' The source code is modified based on code from the blog.
#'
#' @param mat a data/distance matrix.
#' @param diss if `TRUE`, treat `mat` as a distance matrix.
#' @param clsize integer, number of sample within a cluster.
#' @param algo algorithm.
#' @param method method.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set.seed(1234L)
#' x <- rbind(
#'   matrix(rnorm(100, sd = 0.3), ncol = 2),
#'   matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
#' )
#' colnames(x) <- c("x", "y")
#'
#' y1 <- same_size_clustering(x, clsize = 10)
#' y11 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, diss = TRUE)
#'
#' y2 <- same_size_clustering(x, clsize = 10, algo = "hcbottom", method = "ward.D")
#'
#' y3 <- same_size_clustering(x, clsize = 10, algo = "kmvar")
#' y33 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, algo = "kmvar", diss = TRUE)
#' @testexamples
#' expect_length(y1, 100L)
#' expect_length(y11, 100L)
#' expect_length(y2, 100L)
#' expect_length(y3, 100L)
#' expect_length(y33, 100L)
same_size_clustering <- function(mat, diss = FALSE, clsize = NULL,
                                 algo = c("nnit", "hcbottom", "kmvar"),
                                 method = c(
                                   "maxd", "random", "mind", "elki",
                                   "ward.D", "average", "complete", "single"
                                 )) {
  stopifnot(is.numeric(clsize))
  
  algo <- match.arg(algo)
  method <- match.arg(method)
  do.call(algo, args = list(mat = mat, diss = diss, clsize = clsize, method = method))
}

nnit <- function(mat,
                 clsize = NULL,
                 diss = FALSE,
                 method = "maxd") {
  stopifnot(is.logical(diss))
  
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsize <- clsize.rle$lengths
  lab <- rep(NA, nrow(mat))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  cpt <- 1
  while (sum(is.na(lab)) > 0) {
    lab.ii <- which(is.na(lab))
    dmat.m <- dmat[lab.ii, lab.ii]
    ii <- switch(method,
                 maxd = which.max(rowSums(dmat.m)),
                 mind = which.min(rowSums(dmat.m)),
                 random = sample.int(nrow(dmat.m), 1),
                 stop("unsupported method in 'nnit'!")
    )
    lab.m <- rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii, ]), clsize[cpt])] <- cpt
    lab[lab.ii] <- lab.m
    cpt <- cpt + 1
  }
  if (any(is.na(lab))) {
    lab[which(is.na(lab))] <- cpt
  }
  lab
}

kmvar <- function(mat,
                  clsize = NULL,
                  diss = FALSE,
                  method = "maxd") {
  stopifnot(is.logical(diss))
  
  k <- ceiling(nrow(mat) / clsize)
  if (isFALSE(diss)) {
    km.o <- kmeans(mat, k)
    # distance to centers
    centd <- lapply(1:k, function(kk) {
      euc <- t(mat) - km.o$centers[kk, ]
      sqrt(apply(euc, 2, function(x) sum(x^2)))
    })
    centd <- matrix(unlist(centd), ncol = k)
  } else {
    message("PAM algorithm is applied when input distance matrix.")
    pam.o <- cluster::pam(mat, k, diss = TRUE)
    # medoids
    # distance to medoids
    centd <- mat[, pam.o$id.med, drop = FALSE]
  }
  
  labs <- rep(NA, nrow(mat))
  clsizes <- rep(0, k)
  
  ptord <- switch(method,
                  maxd = order(-apply(centd, 1, max)),
                  mind = order(apply(centd, 1, min)),
                  random = sample.int(nrow(mat)),
                  elki = order(apply(centd, 1, min) - apply(centd, 1, max)),
                  stop("unsupported method in 'kmvar'!")
  )
  
  for (ii in ptord) {
    bestcl <- which.max(centd[ii, ])
    labs[ii] <- bestcl
    clsizes[bestcl] <- clsizes[bestcl] + 1
    if (clsizes[bestcl] >= clsize) {
      centd[, bestcl] <- NA
    }
  }
  return(labs)
}

hcbottom <- function(mat,
                     clsize = NULL,
                     diss = FALSE,
                     method = "ward.D") {
  stopifnot(is.logical(diss))
  
  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsizes <- clsize.rle$lengths
  cpt <- 1
  lab <- rep(NA, nrow(mat))
  for (clss in clsizes[-1]) {
    lab.ii <- which(is.na(lab))
    hc.o <- hclust(as.dist(dmat[lab.ii, lab.ii]), method = method)
    clt <- 0
    ct <- length(lab.ii) - clss
    while (max(clt) < clss) {
      cls <- cutree(hc.o, ct)
      clt <- table(cls)
      ct <- ct - 1
    }
    cl.sel <- which(cls == as.numeric(names(clt)[which.max(clt)]))
    lab[lab.ii[head(cl.sel, clss)]] <- cpt
    cpt <- cpt + 1
  }
  lab[is.na(lab)] <- cpt
  lab
}

#' Run MAGICAL
#' 
#' Run MAGICAL
#' 
#' @param loaded_data A list, the output from `prepare_magical_object` or `Data_loading`.
#' @param TAD_file_path The path to the TAD file.
#' @param dc Distance control parameter when TAD file is not specified, Default is 5e5 bps.
#' @param iteration_num The number of iterations for the main estimation step
#' @param Output_file_path The path to write the circuit results
#' @param ... Other parameters for MAGICAL
#' 
#' @return No return value, but will save the inferred circuits and the workspace to files. See the tutorial of MAGICAL.
#' 
#' @export
run_magical_main = function(loaded_data, TAD_file_path, dc = 5e5, iteration_num, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt',...){
  if(exists(TAD_file_path)){
    Candidate_circuits <- Candidate_circuits_construction_with_TAD(loaded_data, TAD_file_path)
  }else{
    Candidate_circuits <- Candidate_circuits_construction_without_TAD(loaded_data, dc)
  }
  Initial_model<-MAGICAL_initialization(loaded_data, Candidate_circuits)
  Circuits_linkage_posterior<-MAGICAL_estimation(loaded_data, Candidate_circuits, Initial_model, iteration_num = 1000)
  MAGICAL_circuits_output(Output_file_path, Candidate_circuits, Circuits_linkage_posterior, ...)
  save(Candidate_circuits, Circuits_linkage_posterior, file = "MAGICAL_results.RData")
}

