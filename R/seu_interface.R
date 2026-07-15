#' Get differential genes and peaks from a Seurat object
#'
#' Run the same contrast on the RNA and ATAC assays of a multiome Seurat
#' object. The returned tables retain the format of [Seurat::FindMarkers()] and
#' can be supplied directly to [Data_loading_from_seu()] as
#' `differential_results`.
#'
#' @param seu A Seurat object containing RNA and ATAC assays.
#' @param ident.1,ident.2 Identities to compare. `ident.2 = NULL` compares
#'   `ident.1` with all other cells, following Seurat's behaviour.
#' @param group.by Optional metadata column used as the identity class for the
#'   comparison.
#' @param subset.ident Optional identity class(es) to subset before applying
#'   `group.by`.
#' @param rna_assay,atac_assay Names of the RNA and ATAC assays.
#' @param p_val_adj Maximum adjusted p-value. Set to `NULL` to skip this filter.
#' @param log2fc Minimum absolute log2 fold change. Set to `NULL` to skip this
#'   filter.
#' @param ... Further arguments passed to [Seurat::FindMarkers()].
#'
#' @return A list with `genes` and `peaks` data frames. Their row names are the
#'   gene and peak identifiers, respectively.
#'
#' @export
Differential_features_from_seu <- function(
    seu,
    ident.1,
    ident.2 = NULL,
    group.by = NULL,
    subset.ident = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    p_val_adj = 0.05,
    log2fc = 0.3,
    ...) {
  .magical_check_seurat_object(seu, rna_assay, atac_assay)

  find_features <- function(assay) {
    markers <- Seurat::FindMarkers(
      object = seu,
      ident.1 = ident.1,
      ident.2 = ident.2,
      group.by = group.by,
      subset.ident = subset.ident,
      assay = assay,
      ...
    )

    fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
    if (length(fc_column) == 0L && !is.null(log2fc)) {
      stop("Seurat did not return an average log-fold-change column.", call. = FALSE)
    }
    keep <- rep(TRUE, nrow(markers))
    if (!is.null(p_val_adj)) {
      if (!"p_val_adj" %in% colnames(markers)) {
        stop("Seurat did not return a p_val_adj column.", call. = FALSE)
      }
      keep <- keep & !is.na(markers$p_val_adj) & markers$p_val_adj <= p_val_adj
    }
    if (!is.null(log2fc)) {
      keep <- keep & !is.na(markers[[fc_column[1L]]]) &
        abs(markers[[fc_column[1L]]]) >= log2fc
    }
    markers[keep, , drop = FALSE]
  }

  list(genes = find_features(rna_assay), peaks = find_features(atac_assay))
}

#' Create MAGICAL input from a Seurat RNA + ATAC co-assay object
#'
#' Extract raw RNA and ATAC counts, feature coordinates, and cell metadata from
#' a multiome Seurat object into the `loaded_data` format consumed by MAGICAL.
#' Candidate genes and peaks may be supplied directly or as the output of
#' [Differential_features_from_seu()]. Motif and RefSeq priors remain external
#' files because they are not generally stored in a Seurat object.
#'
#' @param seu A Seurat object containing RNA and ATAC assays.
#' @param Motif_mapping_file_path A three-column, headerless file containing
#'   ATAC-feature index, motif index, and motif-match value. ATAC indices must
#'   refer to the row order of the ATAC counts layer.
#' @param Motif_name_file_path A two-column, headerless file containing motif
#'   index and motif name.
#' @param Ref_seq_file_path A tab-delimited RefSeq file with columns for
#'   chromosome, strand, start, end, and gene symbol.
#' @param candidate_genes Candidate gene symbols, a one-column data frame, or a
#'   differential-gene table whose row names are gene symbols.
#' @param candidate_peaks Candidate peak coordinates as `chr/start/end` or
#'   `chr/point1/point2` data, a character vector of `chr:start-end` (or
#'   `chr-start-end`) strings, or a differential-peak table whose row names are
#'   peak strings.
#' @param differential_results A list with `genes` and `peaks` components, such
#'   as the return value of [Differential_features_from_seu()]. It is used for
#'   any candidate argument that is `NULL`.
#' @param rna_assay,atac_assay Names of the RNA and ATAC assays.
#' @param cell_type_col,subject_id_col,condition_col Metadata columns used by
#'   MAGICAL. They are copied unchanged into both RNA and ATAC cell metadata.
#' @param count_layer Name of the raw-count layer in each assay.
#'
#' @return A `loaded_data` list suitable for
#'   [Candidate_circuits_construction_with_TAD()] or
#'   [Candidate_circuits_construction_without_TAD()].
#'
#' @export
Data_loading_from_seu <- function(
    seu,
    Motif_mapping_file_path,
    Motif_name_file_path,
    Ref_seq_file_path,
    candidate_genes = NULL,
    candidate_peaks = NULL,
    differential_results = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    cell_type_col = "cell_type",
    subject_id_col = "subject_ID",
    condition_col = "condition",
    count_layer = "counts") {
  .magical_check_seurat_object(seu, rna_assay, atac_assay)

  metadata <- seu[[]]
  metadata_columns <- c(cell_type_col, subject_id_col, condition_col)
  missing_columns <- setdiff(metadata_columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "The Seurat metadata is missing: ", paste(missing_columns, collapse = ", "),
      ". Supply the corresponding *_col arguments.",
      call. = FALSE
    )
  }

  if (!is.null(differential_results)) {
    if (!is.list(differential_results)) {
      stop("`differential_results` must be a list with `genes` and `peaks`.", call. = FALSE)
    }
    if (is.null(candidate_genes)) candidate_genes <- differential_results$genes
    if (is.null(candidate_peaks)) candidate_peaks <- differential_results$peaks
  }
  if (is.null(candidate_genes) || is.null(candidate_peaks)) {
    stop(
      "Supply both candidate arguments or `differential_results`; use ",
      "`Differential_features_from_seu()` to create the latter.",
      call. = FALSE
    )
  }

  rna_counts <- .magical_get_counts(seu[[rna_assay]], count_layer, rna_assay)
  atac_counts <- .magical_get_counts(seu[[atac_assay]], count_layer, atac_assay)
  if (is.null(rownames(rna_counts)) || is.null(rownames(atac_counts))) {
    stop("Both count layers must have feature row names.", call. = FALSE)
  }
  if (is.null(colnames(rna_counts)) || is.null(colnames(atac_counts))) {
    stop("Both count layers must have cell column names.", call. = FALSE)
  }

  scRNA_Genes <- data.frame(
    Gene_index = seq_len(nrow(rna_counts)),
    Gene_symbols = rownames(rna_counts),
    stringsAsFactors = FALSE
  )
  scATAC_Peaks <- .magical_atac_peaks(seu[[atac_assay]], atac_counts)
  Candidate_Genes <- .magical_gene_table(candidate_genes)
  Candidate_Peaks <- .magical_peak_table(candidate_peaks)

  scRNA_cells <- .magical_cell_table(
    metadata, colnames(rna_counts), cell_type_col, subject_id_col, condition_col
  )
  scATAC_cells <- .magical_cell_table(
    metadata, colnames(atac_counts), cell_type_col, subject_id_col, condition_col
  )

  motif_prior <- .magical_read_motif_prior(
    Motif_mapping_file_path, Motif_name_file_path, nrow(atac_counts), rownames(atac_counts)
  )
  Refseq <- utils::read.table(
    Ref_seq_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
    comment.char = "", check.names = FALSE
  )
  if (ncol(Refseq) < 5L) {
    stop("`Ref_seq_file_path` must contain at least five columns.", call. = FALSE)
  }
  Refseq <- Refseq[, seq_len(5L), drop = FALSE]
  colnames(Refseq) <- c("chr", "strand", "start", "end", "Gene_symbols")

  list(
    Common_samples = intersect(unique(scRNA_cells$subject_ID), unique(scATAC_cells$subject_ID)),
    Candidate_Genes = Candidate_Genes,
    Candidate_Peaks = Candidate_Peaks,
    scRNA_Genes = scRNA_Genes,
    scRNA_cells = scRNA_cells,
    scRNA_read_count_matrix = rna_counts,
    scATAC_Peaks = scATAC_Peaks,
    scATAC_cells = scATAC_cells,
    scATAC_read_count_matrix = atac_counts,
    Motifs = motif_prior$motifs,
    TF_Peak_binding_matrix = motif_prior$binding,
    Refseq = Refseq
  )
}

.magical_check_seurat_object <- function(seu, rna_assay, atac_assay) {
  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("`Data_loading_from_seu()` requires the Seurat package.", call. = FALSE)
  }
  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  available_assays <- SeuratObject::Assays(seu)
  missing_assays <- setdiff(c(rna_assay, atac_assay), available_assays)
  if (length(missing_assays) > 0L) {
    stop("The Seurat object is missing assay(s): ", paste(missing_assays, collapse = ", "), call. = FALSE)
  }
}

.magical_get_counts <- function(assay, layer, assay_name) {
  counts <- tryCatch(
    SeuratObject::GetAssayData(assay, layer = layer),
    error = function(e) {
      stop("Could not read the `", layer, "` layer from the ", assay_name,
           " assay: ", conditionMessage(e), call. = FALSE)
    }
  )
  methods::as(counts, "dgCMatrix")
}

.magical_cell_table <- function(metadata, cells, cell_type_col, subject_id_col, condition_col) {
  absent_cells <- setdiff(cells, rownames(metadata))
  if (length(absent_cells) > 0L) {
    stop("Metadata is missing cells present in an assay.", call. = FALSE)
  }
  metadata <- metadata[cells, , drop = FALSE]
  data.frame(
    cell_index = seq_along(cells),
    cell_barcode = cells,
    cell_type = as.character(metadata[[cell_type_col]]),
    subject_ID = as.character(metadata[[subject_id_col]]),
    condition = as.character(metadata[[condition_col]]),
    stringsAsFactors = FALSE
  )
}

.magical_gene_table <- function(genes) {
  if (is.data.frame(genes) || is.matrix(genes)) {
    genes <- if (ncol(genes) == 1L) genes[, 1L] else rownames(genes)
  }
  if (is.null(genes) || is.null(genes <- as.character(genes))) {
    stop("Candidate genes must have gene symbols.", call. = FALSE)
  }
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0L) stop("No candidate genes were supplied.", call. = FALSE)
  data.frame(Gene_symbols = genes, stringsAsFactors = FALSE)
}

.magical_peak_table <- function(peaks) {
  if (inherits(peaks, "GRanges")) {
    peaks <- data.frame(
      chr = as.character(GenomicRanges::seqnames(peaks)),
      point1 = GenomicRanges::start(peaks),
      point2 = GenomicRanges::end(peaks),
      stringsAsFactors = FALSE
    )
  }
  if (is.data.frame(peaks) || is.matrix(peaks)) {
    peak_names <- colnames(peaks)
    if (all(c("chr", "point1", "point2") %in% peak_names)) {
      peaks <- peaks[, c("chr", "point1", "point2"), drop = FALSE]
    } else if (all(c("chr", "start", "end") %in% peak_names)) {
      peaks <- peaks[, c("chr", "start", "end"), drop = FALSE]
    } else if (!is.null(rownames(peaks))) {
      peaks <- rownames(peaks)
    } else {
      stop("Candidate peaks need coordinate columns or peak-string row names.", call. = FALSE)
    }
  }
  if (is.character(peaks)) {
    matched <- regexec("^(.+?)(?::|-)([0-9]+)-([0-9]+)$", peaks)
    parts <- regmatches(peaks, matched)
    if (any(lengths(parts) != 4L)) {
      stop("Peak strings must use `chr:start-end` or `chr-start-end`.", call. = FALSE)
    }
    peaks <- data.frame(
      chr = vapply(parts, `[`, character(1L), 2L),
      point1 = as.numeric(vapply(parts, `[`, character(1L), 3L)),
      point2 = as.numeric(vapply(parts, `[`, character(1L), 4L)),
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(peaks) && !is.matrix(peaks)) {
    stop("Candidate peaks must be coordinates, peak strings, or a GRanges object.", call. = FALSE)
  }
  peaks <- as.data.frame(peaks, stringsAsFactors = FALSE)
  if (ncol(peaks) != 3L) stop("Candidate peaks must have three coordinate columns.", call. = FALSE)
  colnames(peaks) <- c("chr", "point1", "point2")
  peaks$chr <- as.character(peaks$chr)
  peaks$point1 <- suppressWarnings(as.numeric(peaks$point1))
  peaks$point2 <- suppressWarnings(as.numeric(peaks$point2))
  if (anyNA(peaks) || any(peaks$point1 > peaks$point2)) {
    stop("Candidate peak coordinates must be complete and have start <= end.", call. = FALSE)
  }
  unique(peaks)
}

.magical_atac_peaks <- function(atac_assay, atac_counts) {
  ranges <- tryCatch(Signac::granges(atac_assay), error = function(e) NULL)
  if (!is.null(ranges) && length(ranges) == nrow(atac_counts)) {
    return(data.frame(
      Peak_index = seq_len(nrow(atac_counts)),
      chr = as.character(GenomicRanges::seqnames(ranges)),
      point1 = GenomicRanges::start(ranges),
      point2 = GenomicRanges::end(ranges),
      stringsAsFactors = FALSE
    ))
  }
  parsed <- .magical_peak_table(rownames(atac_counts))
  data.frame(Peak_index = seq_len(nrow(atac_counts)), parsed, stringsAsFactors = FALSE)
}

.magical_read_motif_prior <- function(mapping_path, motif_path, n_peaks, peak_names) {
  motifs <- utils::read.table(
    motif_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE
  )
  mapping <- utils::read.table(
    mapping_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE
  )
  if (ncol(motifs) < 2L || ncol(mapping) < 3L) {
    stop("Motif names and motif mapping files must have at least two and three columns, respectively.", call. = FALSE)
  }
  motifs <- motifs[, 1:2, drop = FALSE]
  mapping <- mapping[, 1:3, drop = FALSE]
  motif_ids <- as.character(motifs[[1L]])
  motif_positions <- match(as.character(mapping[[2L]]), motif_ids)
  if (anyNA(motif_positions)) {
    stop("Motif indices in the mapping file do not match the motif-name file.", call. = FALSE)
  }
  peak_positions <- suppressWarnings(as.integer(mapping[[1L]]))
  values <- suppressWarnings(as.numeric(mapping[[3L]]))
  if (anyNA(peak_positions) || anyNA(values) || any(peak_positions < 1L | peak_positions > n_peaks)) {
    stop("The motif mapping has invalid ATAC-feature indices or match values.", call. = FALSE)
  }
  list(
    motifs = data.frame(motif_index = seq_len(nrow(motifs)), name = as.character(motifs[[2L]]), stringsAsFactors = FALSE),
    binding = Matrix::sparseMatrix(
      i = peak_positions,
      j = motif_positions,
      x = values,
      dims = c(n_peaks, nrow(motifs)),
      dimnames = list(peak_names, as.character(motifs[[2L]]))
    )
  )
}
