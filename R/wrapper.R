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
#' @param to_file Whether to write the MAGICAl inputs to files under the folder "input files".
#' @param genome The genome for searching TF binding motifs. Default is `"hg38"`.
#' @param meta_spot_opt To run MAGICAL at meta-spot level or not. Default is `F`.
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' @param cl_method A choice from "mclust" and "louvain". The clustering method.
#' @param mclust.num If `method = "mclust`, this is the number of clusters.
#' @param ld.resolution If `method = "louvain"`, this is the resolution parameter.
#' @param random.seed Random seed.
#' @param TAD_file_path The path to the TAD file.
#' @param dc Distance control parameter when TAD file is not specified, Default is 5e5 bps.
#' @param iteration_num The number of iterations for the main estimation step
#' @param Output_file_path The path to write the circuit results
#' @param ... Other parameters
#' 
#' @return If `magical = F`: a list of filtered data.frame in the form of the result of `Seurat::FindMarkers()`. If `magical = T`: a list of 2 data.frames: 1 is the one above, the other is the output from `MAGICAL_circuits_output()`.
wrapper_main = function(RNA_counts, ATAC_counts, niche_label, meta_add=NULL, pb, contrast = c("niche","condition"), niche1=NULL, niche2 = NULL, condition=NULL, condition1 = NULL, condition2 = NULL, p_thre = 0.05, log2fc_thre = 0.3, magical, to_file = F, Ref_seq_file_path, genome = "hg38", meta_spot_opt = F, feature,cl_method, mclust.num, ld.resolution, random.seed, TAD_file_path, dc = 5e5, iteration_num = 250, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt', ...){
  ## step 1: differential analysis
  cat("Performing differential analysis ... \n")
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
    if(nsample < 10){
      if(!exists(meta_spot_opt)){
        meta_spot_opt = T
        print("There are less than 10 samples. Meta-spot level MAGICAL analysis will be applied.")
      }
      else{
        if(meta_spot_opt == F){
          warning("There are less than 10 samples. We suggest using meta-spot level MAGICAL analysis.")
        }
      }
    }else{
      if(!exists(meta_spot_opt)){
        meta_spot_opt = F
        print("There are more than 10 samples. Pseudo-bulk level differential analysis will be applied. You can change to meta-spot level by setting `meta_spot = T`.")
      }
    }
    
    ## In this version, only `swk` clustering is retained. It will use 6 different sets of clustering settings and run 6 MAGICAL separately. Accordingly, `iteration_num` is downgraded.
    
    cl_settings = data.frame(cl_method = rep(c("mclust","louvain"),each = 3), param = c(5,10,20,0.5,2,4))
    circuits = list()
    
    for(i in 1:3){
      ## step 2: prepare MAGICAL input
      cat("\n Preparing inputs for MAGICAL ...\n")
      loaded_data = prepare_magical_object(differentials[["deg"]], differentials[["das"]], RNA_counts, ATAC_counts, niche_label, meta_add, Ref_seq_file_path, meta_spot_opt, contrast, niche1, niche2, condition, condition1, condition2, feature, cl_method = "mclust", mclust.num = cl_settings[i,2], random.seed, ...)
      
      ## step 3: run MAGICAL
      cat("\n Running MAGICAL ... \n")
      run_magical_main(loaded_data, TAD_file_path, dc, iteration_num, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt',...)
    }
    for(i in 4:6){
      cat("\n Preparing inputs for MAGICAL ...\n")
      loaded_data = prepare_magical_object(differentials[["deg"]], differentials[["das"]], RNA_counts, ATAC_counts, niche_label, meta_add, Ref_seq_file_path, meta_spot_opt, contrast, niche1, niche2, condition, condition1, condition2, feature, cl_method = "louvain", ld.resolution = cl_settings[i,2], random.seed, ...)
      
      cat("\n Running MAGICAL ... \n")
      circuits[[i]] = run_magical_main(loaded_data, TAD_file_path, dc, iteration_num, Output_file_path = 'MAGICAL_selected_regulatory_circuits.txt',...)
    }
    
    ## find consensus circuits
    consensus_circuits = Reduce(function(x, y) inner_join(x, y, by = c("Gene_symbol", "Gene_chr", "Gene_TSS", "Peak_chr", "Peak_start", "Peak_end")), circuits)
  }
  return(list(differentials,consensus_circuits))
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
    cat("Building pseudo-bulks ... \n")
    
    spots_niche1id = which(niche_label == niche1)
    if(contrast =="niche"){
      RNA_niche1 = RNA_counts[,spots_niche1id]
      ATAC_niche1 = ATAC_counts[,spots_niche1id]
      meta_niche1 = meta_add[spots_niche1id,"sample"]
      sample_df1 = data.frame(spot = rownames(meta_add)[spots_niche1id], sample = meta_niche1)
      if(is.null(niche2)){
        niche2 = paste0("non_",niche1)
        RNA_niche2 = RNA_counts[,-spots_niche1id]
        ATAC_niche2 = ATAC_counts[,-spots_niche1id]
        meta_niche2 = meta_add[-spots_niche1id,"sample"]
        sample_df2 = data.frame(spot = rownames(meta_add)[-spots_niche1id], sample = meta_niche2)
      }else{
        spots_niche2id = which(niche_label == niche2)
        RNA_niche2 = RNA_counts[,spots_niche2id]
        ATAC_niche2 = ATAC_counts[,spots_niche2id]
        meta_niche2 = meta_add[spots_niche2id,"sample"]
        sample_df2 = data.frame(spot = rownames(meta_add)[spots_niche2id], sample = meta_niche2)
      }
      
      # aggregate to pb
      RNA_agg1 = aggregate_to_pb(sample_df1, RNA_niche1)
      ATAC_agg1 = aggregate_to_pb(sample_df1,ATAC_niche1)
      
      RNA_agg2 = aggregate_to_pb(sample_df2, RNA_niche2)
      ATAC_agg2 = aggregate_to_pb(sample_df2,ATAC_niche2)
    }else if(contrast == "condition"){
      RNA_counts_totake = RNA_counts[,spots_niche1id]
      ATAC_counts_totake = ATAC_counts[,spots_niche1id]
      niche_label_totake = niche_label[spots_niche1id]
      meta_add_totake = meta_add[spots_niche1id,]
      
      spots_condition1id = which(meta_add_totake[[condition]] == condition1)
      
      RNA_condition1 = RNA_counts_totake[,spots_condition1id]
      ATAC_condition1 = ATAC_counts_totake[,spots_condition1id]
      meta_condition1 = meta_add[spots_condition1id,"sample"]
      sample_df1 = data.frame(spot = rownames(meta_add_totake)[spots_condition1id], sample = meta_condition1)
      if(is.null(condition2)){
        condition2 = paste0("non_",condition1)
        RNA_condition2 = RNA_counts_totake[,-spots_condition1id]
        ATAC_condition2 = ATAC_counts_totake[,-spots_condition1id]
        meta_condition2 = meta_add_totake[-spots_condition1id,"sample"]
        sample_df2 = data.frame(spot = rownames(meta_add_totake)[-spots_condition1id], sample = meta_condition2)
      }else{
        spots_condition2id = which(meta_add[[condition]] == condition2)
        RNA_condition2 = RNA_counts_totake[,spots_condition2id]
        ATAC_condition2 = ATAC_counts_totake[,spots_condition2id]
        meta_condition2 = meta_add_totake[spots_condition2id,"sample"]
        sample_df2 = data.frame(spot = rownames(meta_add_totake)[spots_condition2id], sample = meta_condition2)
      }
      
      # aggregate to pb
      RNA_agg1 = aggregate_to_pb(sample_df1, RNA_condition1)
      ATAC_agg1 = aggregate_to_pb(sample_df1,ATAC_condition1)
      
      RNA_agg2 = aggregate_to_pb(sample_df2, RNA_condition2)
      ATAC_agg2 = aggregate_to_pb(sample_df2,ATAC_condition2)
    }
    
    colnames(RNA_agg2) = paste0(colnames(RNA_agg2),"_2")
    colnames(ATAC_agg2) = paste0(colnames(ATAC_agg2),"_2")
    RNA = cbind(RNA_agg1, RNA_agg2)
    ATAC = cbind(ATAC_agg1, ATAC_agg2)
    
    # new metadata
    metadata = data.frame(
      spot = colnames(RNA),
      sample = c(unique(meta_niche1),unique(meta_niche2)),
      niche_label = c(colnames(RNA_agg1),colnames(RNA_agg2))
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
    cells1 = obj$cell[which(obj$niche_label == niche1)]
    if(is.null(niche2)){
      cells2 = obj$cell[which(obj$niche_label != niche1)]
    }else{cells2 = obj$cell[which(obj$niche_label == niche2)]}
  }else if(contrast == "condition"){
    cells1 = obj$cell[which(obj$niche_label == niche1 & obj[[conditionname]] == condition1)]
    if(is.null(condition2)){
      cells2 = obj$cell[which(obj$niche_label == niche1 & obj[[conditionname]] != condition1)]
    }else{obj$cell[which(obj$niche_label == niche1 & obj[[conditionname]] == condition2)]}
  }
  
  cat("Selecting differentially expressed genes ...\n")
  diff_genes = Seurat::FindMarkers(obj@assays$RNA, cells.1 = cells1, cells.2 = cells2, ...)
  deg = diff_genes[which(diff_genes$p_val_adj<p_thre & abs(diff_genes$avg_log2FC)>log2fc_thre),]
  paste0("We get ", dim(deg)[1], " genes.")
  
  cat("Selecting differentially associated sites ... \n")
  diff_peaks = Seurat::FindMarkers(obj@assays$ATAC, cells.1 = cells1, cells.2 = cells2, ...)
  das = diff_peaks[which(diff_peaks$p_val_adj<p_thre & abs(diff_peaks$avg_log2FC)>log2fc_thre),]
  paste0("We get ", dim(das)[1], " peaks.")
  
  return(list(genes = deg, peaks = das))
}

aggregate_to_pb = function(meta, mtx){
  facts = unique(meta$sample)
  agg_mtx = matrix(nrow = nrow(mtx), ncol = length(facts))
  colnames(agg_mtx) = facts
  rownames(agg_mtx) = rownames(mtx)
  for(fact in facts){
    totake = meta[which(meta$sample %in% fact),"spot"]
    if(length(totake)==1){
      agg_mtx[,fact] = mtx[,totake]
    }else{
      toadd = mtx[,totake]
      agg_mtx[,fact] = rowSums(toadd)
    }
  }
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
#' @param to_file Whether to write the MAGICAl inputs to files under the folder "input files".
#' @param Ref_seq_file_path Path to the Refseq file for transcription starting site extraction
#' @param genome The genome for searching TF binding motifs. Default is `"hg38"`.
#' @param meta_spot_opt To run MAGICAL at meta-spot level or not. Default is `F`.
#' @param contrast Niche- or condition-specific contrast
#' @param niche1 (Required for all cases) For `contrast = "niche"`, this should be one niche you want to contrast. For `contrast = "condition"`, this should be the niche in which you want to contrast the conditions
#' @param niche2 The other niche you want to contrast if `contrast = "niche"`. Not specified is to contrast `niche1` with all other niches
#' @param conditionname The condition you want to make the contrast if `contrast = "condition"`
#' @param condition1 One condition you want to contrast if `contrast = "condition"`
#' @param condition2 The other condition you want to contrast if `contrast = "condition"`. Not specified is to contrast `condition1` with all other conditions.
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' @param cl_method A choice from "mclust" and "louvain". The clustering method for `swk()`.
#' @param mclust.num If `method = "mclust`, this is the number of clusters.
#' @param ld.resolution If `method = "louvain"`, this is the resolution parameter.
#' @param random.seed Random seed.
#' @param ... Other parameters for `same_size_clustering()` or `swk()`
#' 
#' @import Matrix
#' @import Seurat
#' @import Signac
#' @import TFBSTools
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import chromVARmotifs
#'
#' @export
prepare_magical_object = function(deg, das, RNA_counts, ATAC_counts, niche_label, meta_add, to_file = F, Ref_seq_file_path, genome = "hg38", meta_spot_opt = F, contrast = c("niche","condition"), niche1=NULL, niche2 = NULL, condition=NULL, condition1 = NULL, condition2 = NULL, feature, cl_method, mclust.num, ld.resolution, random.seed, ...){
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
  
  candidate_peaks_char = rownames(das)
  candidate_peaks = do.call(rbind, strsplit(candidate_peaks_char, "-"))
  candidate_peaks = data.frame(chr = candidate_peaks[,1], point1 = as.numeric(candidate_peaks[,2], point2 = as.numeric(candidate_peaks[,3])))
  
  RNA_count_mtx = as(RNA_counts, "TsparseMatrix")
  
  RNA_genes = data.frame(Gene_index = 1:nrow(RNA_counts), Gene_symbols = rownames(RNA_counts))
  
  ATAC_count_mtx = as(ATAC_counts, "TsparseMatrix")
  
  ATAC_peaks_char = rownames(ATAC_counts)
  ATAC_peaks = Signac::StringToGRanges(ATAC_peaks_char, sep = c(":", "-"))
  
  Refseq = read.table(Ref_seq_file_path, header = TRUE, sep = "\t")
  colnames(Refseq) = c("chr", "strand", "start", "end", "Gene_symbols")
  
  # get motifs with chromVARmotifs
  cat("\n Getting motfis. This step make take some time ... \n")
  chromatinassay <- Signac::CreateChromatinAssay(counts = ATAC_counts, genome = genome)
  object <- Seurat::CreateSeuratObject(counts = chromatinassay)
  library("BSgenome.Hsapiens.UCSC.hg38") ###
  object <- Signac::AddMotifs(object = object, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = human_pwms_v2) ###this line need to be updated (e.g., if genome == "hg38", ...)
  Peak_motif_mapping <- as(object@assays$RNA@motifs@data * 1, "TsparseMatrix")
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
  colnames(motifs) <- c("motif_index", "name")
  
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
    cat("Building meta-spots ... \n")
    if(contrast == "niche"){
      clusters = niche_label
      if(is.null(niche2)){cluster[which(cluster!=niche1)] = paste0("non_",niche1)}
    }else if(contrast == "condition"){
      clusters = meta_add[[condition]]
      if(is.null(condition2)){cluster[which(cluster!=condition1)] = paste0("non_",condition1)}
    }
    
    new_idents = metaspot_swk(clusters, feature, cl_method, mclust.num, ld.resolution, random.seed, ...)
    if(contrast == "niche"){
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche_label, subject_ID = new_idents, condition = "1")
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche_label, subject_ID = new_idents, condition = "1")
    }else if(contrast == "condition"){
      RNA_cells = data.frame(cell_index = 1:ncol(RNA_counts), cell_barcode = colnames(RNA_counts), cell_type = niche1, subject_ID = new_idents, condition = meta_add[[condition]])
      ATAC_cells = data.frame(cell_index = 1:ncol(ATAC_counts), cell_barcode = colnames(ATAC_counts), cell_type = niche1, subject_ID = new_idents, condition = meta_add[[condition]])
    }
  }
  
  Common_samples <- intersect(RNA_cells$subject_ID, ATAC_cells$subject_ID)
  
  if(to_file == T){
    cat("Writing files ...\n")
    if(!dir.exists("input files"))dir.create("input_files")
    write.table(candidate_genes, file = "input files/Cell type candidate genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(candidate_peaks, file = "input files/Cell type candidate peaks.txt", quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(summary(RNA_count_mtx), file = "input files/Cell type scRNA read count.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(RNA_genes, file = "input files/scRNA genes.txt", quote = F, col.names = F, sep = "\t")
    write.table(RNA_cells, file = "input files/Cell type scRNA cell meta.txt", quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(summary(ATAC_count_mtx), file = "/data/home/zz5708/Projects/KPMP_V2/MATLAB/input files/Cell type scATAC read count.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(ATAC_peaks, file = "input files/scATAC peaks.txt", quote = F, col.names = F, sep = "\t")
    write.table(ATAC_cells, file = "input files/Cell type scATAC cell meta.txt", quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(motif_mapping, file = "input files/Motif mapping prior.txt", quote = F, row.names = F, col.names = F, sep = "\t")
    write.table(motifs, file = "input files/Motifs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  }
  
  loaded_data = list(
    "Common_samples" = Common_samples,
    "Candidate_Genes" = data.frame(candidate_genes),
    "Candidate_Peaks" = data.frame(candidate_peaks),
    "scRNA_Genes" = data.frame(RNA_genes),
    "scRNA_cells" = data.frame(RNA_cells),
    "scRNA_read_count_matrix" = RNA_count_mtx,
    "scATAC_Peaks" = data.frame(ATAC_peaks),
    "scATAC_cells" = data.frame(ATAC_cells),
    "scATAC_read_count_matrix" = ATAC_count_mtx,
    "Motifs" = data.frame(motifs),
    "TF_Peak_binding_matrix" = motif_mapping,
    "Refseq" = data.frame(Refseq)
  )
  
  return(loaded_data)
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
#' @return The return value of `MAGICAL_circuits_output()`.
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
  res = MAGICAL_circuits_output(Output_file_path, Candidate_circuits, Circuits_linkage_posterior, ...)
  # save(Candidate_circuits, Circuits_linkage_posterior, file = "MAGICAL_results.RData")
  return(res)
}

#' Functions for metaspot construction
#' 
#' Using `DAVINCI::swk()`
#' 
#' @param clusters The labels of the spots
#' @param feature The features to be clustered. Can be LVs from DAVINCI, or spatial coordinates.
#' @param cl_method A choice from "mclust" and "louvain". The clustering method.
#' @param mclust.num If `method = "mclust`, this is the number of clusters.
#' @param ld.resolution If `method = "louvain"`, this is the resolution parameter.
#' @param random.seed Random seed.
#' @param ... Other parameters required by `DAVINCI::swk()`.
#' 
#' @export
metaspot_swk = function(clusters, feature, cl_method = c("mclust", "louvain"), mclust.num = NULL, ld.resolution = NULL, random.seed = 1, ...) {
  results <- data.frame(index = integer(), label = character(), group_index = integer())
  labels <- unique(clusters)
  
  for (label in labels) {
    totake <- which(clusters == label)
    if(length(totake)==0){
      next()
    }
    
    subset_clusters <- cbind.data.frame(1:length(totake), clusters[totake])
    
    subset_features <- feature[subset_clusters[,1], ]
    sink("/dev/null") # to suppress messages from `print()`
    result <-  suppressMessages(swk(subset_features, cl_method, mclust.num = mclust.num, ld.resolution = ld.resolution, random.seed, ...))
    subset_clusters$group_index <- paste0(label, "-", result)
    sink()
    
    results <- rbind(results, subset_clusters)
  }
  return(results)
}