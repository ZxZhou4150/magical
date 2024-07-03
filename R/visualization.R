#' Plot circuit(s)
#'
#' This function plots one circuit by its index.
#'
#' @seealso [plot_circuits_with_gene()] and [plot_circuits_with_peak()]
#'
#' @param data The data frame written by [MAGICAL_circuits_output()].
#' @param idx The index of the circuit to be plotted.
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `TRUE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#' @param peak_track Whether to plot peak tracks around the region or not. Default is `TRUE`.
#' @param peaks The data frame containing the peaks. It should contain 3 columns: the chromosome, start, and end of the peaks. Required if `peak_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' # data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #
#' # library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #
#' # peaks = read.table('Demo input files/scATAC peaks.txt', header = F, sep = "\t")
#' #
#' # plot_circuits_with_idx(data, 1, gene_track = T, TxDb.Hsapiens.UCSC.hg38.knownGene, peak_track = T, peaks)
#'
#' @export

plot_circuits_with_idx <- function(data, idx, gene_track = T, TxDb, peak_track = T, peaks) {
  if (length(idx) != 1) {
    stop("idx should be a single number")
  }

  data_toplot <- data[idx, ]

  ## finding region
  chr <- data_toplot$Gene_chr
  maxsite <- max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite <- min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit <- floor(log10(maxsite - minsite))
  start <- minsite - minsite %% 10^length_digit - 10^length_digit / 2
  end <- maxsite + (-maxsite %% 10^length_digit) + 10^length_digit / 2

  region <- regioneR::toGRanges(paste0(chr, ":", start, "-", end))
  pp <- karyoploteR::getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 2
  kp <- karyoploteR::plotKaryotype(cytobands = GenomicRanges::GRanges(), zoom = region, plot.params = pp, plot.type = 2, cex = 1)

  ## gene track
  if (gene_track) {
    if (missing(TxDb)) {
      stop("TxDb is required. Please install the corresponding TxDb package")
    }
    genes.data <- karyoploteR::makeGenesDataFromTxDb(TxDb, karyoplot = kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data = genes.data, data.panel = 2, r0 = 0.7, r1 = 0.85, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
    karyoploteR::kpAddLabels(kp, labels = "Gene track", data.panel = 2, r0 = 0.7, r1 = 0.85, cex = 0.6, col = "black")
  }

  ## peak track
  if (peak_track) {
    if (missing(peaks)) {
      stop("peaks are required")
    }
    peaks_range <- regioneR::toGRanges(peaks)
    # karyoploteR::kpPlotRegions(kp, data = peaks_range, data.panel = 2, r0 = 0.9, r1 = 1, col = peaks_range$itemRgb)
    karyoploteR::kpPlotDensity(kp, data = peaks_range, data.panel = 2, r0 = 0.4, r1 = 0.1, window.size = 10^(length_digit - 1.1))
    karyoploteR::kpAddLabels(kp, labels = "Chromatin \n activity", data.panel = 2, r0 = 0.4, r1 = 0.1, cex = 0.6, col = "black")
  }

  ## finding start and end

  starts <- regioneR::toGRanges(data_toplot[, c("Gene_chr", "Gene_TSS", "Gene_TSS")])
  ends <- regioneR::toGRanges(data_toplot[, c("Peak_chr", "Peak_start", "Peak_end")])

  ## plotting
  karyoploteR::kpPlotRegions(kp, data = starts, r0 = 0.2, r1 = 0.3, col = "blue")
  karyoploteR::kpPlotRegions(kp, data = ends, r0 = 0.2, r1 = 0.3, col = "red")
  karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.3, r1 = 0.5, col = "grey")
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit - 1), digit = 6 - length_digit, add.units = TRUE, cex = 0.6, tick.len = 3)

  karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[, "Gene_TSS"], y = 0.05, col = "black", cex = 0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.1, y1 = 0.18, length = 0.05)

  karyoploteR::kpText(kp, labels = paste0("peak ", chr, ":", data_toplot$Peak_start, "-", data_toplot$Peak_end), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y = 0.65, col = "black", cex = 0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, x1 = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y0 = 0.5, y1 = 0.4, length = 0.05)

  # labels <- paste0("peak ", data_toplot$Peak_chr, ":", data_toplot$Peak_start, "-", data_toplot$Peak_end)
  # peak_marker <- data.frame(chr = data_toplot$Peak_chr, pos = (data_toplot$Peak_start + data_toplot$Peak_end) / 2,labels = labels)
  # karyoploteR::kpPlotMarkers(kp,chr = peak_marker$chr, x = peak_marker$pos, y = 0.65, labels = peak_marker$labels,text.orientation = "horizontal",adjust.label.position = T, line.color = "white", cex = 0.6)

  ### top TF(s)
  TF_binding_prob <- data_toplot[, "TFs.binding.prob."]
  TFs <- gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
  TFs <- TFs[TFs != ""]
  karyoploteR::kpText(kp, labels = paste0("top binding TF(s): ", paste(head(TFs, 3), collapse = " ")), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y = 0.55, col = "black", cex = 0.6)

  ### title
  bb <- karyoploteR::getMainTitleBoundingBox(kp)
  x <- (bb$x0 + bb$x1) / 2
  y <- 0.9

  graphics::text(x = x, y = y, labels = paste0("Circuit #", idx), cex = 1, font = 2)
}

#' Plot circuit(s)
#'
#' This function plots all circuits containing one certain gene.
#'
#' @seealso [plot_circuits_with_idx()] and [plot_circuits_with_peak()]
#'
#' @param data The data frame written by [MAGICAL_circuits_output()].
#' @param gene The gene symbol.
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `TRUE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#' @param peak_track Whether to plot peak tracks around the region or not. Default is `TRUE`.
#' @param peaks The data frame containing the peaks. It should contain 3 columns: the chromosome, start, and end of the peaks. Required if `peak_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' # data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #
#' # library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #
#' # peaks = read.table('Demo input files/scATAC peaks.txt', header = F, sep = "\t")
#' #
#' # plot_circuits_with_gene(data, gene = "ACOT11", gene_track = T, TxDb.Hsapiens.UCSC.hg38.knownGene, peak_track = T, peaks)
#'
#' @export

plot_circuits_with_gene <- function(data, gene, gene_track = T, TxDb, peak_track = T, peaks) {
  if (length(gene) != 1) {
    stop("gene should be a single gene symbol")
  }
  idx <- which(data$Gene_symbol == gene)
  if (length(idx) == 0) {
    stop("gene not found")
  }

  data_toplot <- data[idx, ]

  ## finding region
  chr <- data_toplot[1, "Gene_chr"]
  maxsite <- max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite <- min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit <- floor(log10(maxsite - minsite))
  start <- minsite - minsite %% 10^length_digit - 10^length_digit / 2
  end <- maxsite + (-maxsite %% 10^length_digit) + 10^length_digit / 2

  region <- regioneR::toGRanges(paste0(chr, ":", start, "-", end))
  pp <- karyoploteR::getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 2
  kp <- karyoploteR::plotKaryotype(cytobands = GenomicRanges::GRanges(), zoom = region, plot.params = pp, plot.type = 2, cex = 1)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit - 1), digit = 6 - length_digit, add.units = TRUE, cex = 0.6, tick.len = 3)

  ## gene track
  if (gene_track) {
    if (missing(TxDb)) {
      stop("TxDb is required. Please install the corresponding TxDb package")
    }
    genes.data <- karyoploteR::makeGenesDataFromTxDb(TxDb, karyoplot = kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data = genes.data, data.panel = 2, r0 = 0.7, r1 = 0.85, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
    karyoploteR::kpAddLabels(kp, labels = "Gene track", data.panel = 2, r0 = 0.7, r1 = 0.85, cex = 0.6, col = "black")
  }

  ## peak track
  if (peak_track) {
    if (missing(peaks)) {
      stop("peaks are required")
    }
    peaks_range <- regioneR::toGRanges(peaks)
    # karyoploteR::kpPlotRegions(kp, data = peaks_range, data.panel = 2, r0 = 0.9, r1 = 1, col = peaks_range$itemRgb)
    karyoploteR::kpPlotDensity(kp, data = peaks_range, data.panel = 2, r0 = 0.4, r1 = 0.1, window.size = 10^(length_digit - 1.1))
    karyoploteR::kpAddLabels(kp, labels = "Chromatin \n activity", data.panel = 2, r0 = 0.4, r1 = 0.1, cex = 0.6, col = "black")
  }

  ## plot gene
  starts <- regioneR::toGRanges(data_toplot[1, c("Gene_chr", "Gene_TSS", "Gene_TSS")])
  karyoploteR::kpPlotRegions(kp, data = starts, r0 = 0.2, r1 = 0.3, col = "blue")
  karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[, "Gene_TSS"], y = 0.05, col = "black", cex = 0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.1, y1 = 0.18, length = 0.05)

  ## plot peaks
  labels <- paste0("peak ", data_toplot$Peak_chr, ":", data_toplot$Peak_start, "-", data_toplot$Peak_end)

  for (i in 1:nrow(data_toplot)) {
    ends <- regioneR::toGRanges(data_toplot[i, c("Peak_chr", "Peak_start", "Peak_end")])
    karyoploteR::kpPlotRegions(kp, data = ends, r0 = 0.2, r1 = 0.3, col = "red")
    karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.3, r1 = 0.5, col = "grey")
    # karyoploteR::kpText(kp, labels = paste0("peak #", i, " ", chr, ":", data_toplot[i, "Peak_start"], "-", data_toplot[i, "Peak_end"]), chr = chr, x = (data_toplot[i, "Peak_start"] + data_toplot[i, "Peak_end"]) / 2, y = 0.65 - 0.05 * (i - 1), col = "black", cex = 0.6)
    karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot[i, "Peak_start"] + data_toplot[i, "Peak_end"]) / 2, x1 = (data_toplot[i, "Peak_start"] + data_toplot[i, "Peak_end"]) / 2, y0 = 0.45, y1 = 0.4, length = 0.05)


    TF_binding_prob <- data_toplot[i, "TFs.binding.prob."]
    TFs <- gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
    TFs <- TFs[TFs != ""]

    labels[i] <- paste(labels[i], paste0("top binding TF(s): ", paste(head(TFs, 3), collapse = " ")), sep = "\n")
  }

  ### top TF(s)

  peak_marker <- data.frame(chr = data_toplot$Peak_chr, pos = (data_toplot$Peak_start + data_toplot$Peak_end) / 2,labels = labels)
  karyoploteR::kpPlotMarkers(kp,chr = peak_marker$chr, x = peak_marker$pos, y = 0.6, labels = peak_marker$labels,text.orientation = "horizontal",adjust.label.position = T, line.color = "black", cex = 0.6, r0 = 0.4, r1 = 0.6)

  ### title
  bb <- karyoploteR::getMainTitleBoundingBox(kp)
  x <- (bb$x0 + bb$x1) / 2
  y <- 0.9

  graphics::text(x = x, y = y, labels = paste0("Circuits containing gene ", gene), cex = 1, font = 2)
}

#' Plot circuit(s)
#'
#' This function plots all circuits containing one certain peak.
#'
#' @seealso [plot_circuits_with_idx()] and [plot_circuits_with_gene()]
#'
#' @param data The data frame written by [MAGICAL_circuits_output()].
#' @param peak_chr (Character) the chromosome of the peak.
#' @param peak_start The start position of the peak.
#' @param peak_end The end position of the peak.
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `TRUE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#' @param peak_track Whether to plot peak tracks around the region or not. Default is `TRUE`.
#' @param peaks The data frame containing the peaks. It should contain 3 columns: the chromosome, start, and end of the peaks. Required if `peak_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' # data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #
#' # library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #
#' # peaks = read.table('Demo input files/scATAC peaks.txt', header = F, sep = "\t")
#' #
#' # plot_circuits_with_peak(data, peak_chr = "chr4", peak_start = 184817550, peak_end = 184818412, gene_track = T, TxDb.Hsapiens.UCSC.hg38.knownGene, peak_track = T, peaks)
#'
#' @export

plot_circuits_with_peak <- function(data, peak_chr, peak_start, peak_end, gene_track = T, TxDb, peak_track = T, peaks) {
  idx <- which(data$Peak_chr == peak_chr & data$Peak_start == peak_start & data$Peak_end == peak_end)
  if (length(idx) == 0) {
    stop("peak not found")
  }

  data_toplot <- data[idx, ]

  ## finding region
  chr <- peak_chr
  maxsite <- max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite <- min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit <- floor(log10(maxsite - minsite))
  start <- minsite - minsite %% 10^length_digit - 10^length_digit / 2
  end <- maxsite + (-maxsite %% 10^length_digit) + 10^length_digit / 2

  region <- regioneR::toGRanges(paste0(chr, ":", start, "-", end))
  pp <- karyoploteR::getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <- 2
  kp <- karyoploteR::plotKaryotype(cytobands = GenomicRanges::GRanges(), zoom = region, plot.params = pp, plot.type = 2, cex = 1)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit - 1), digit = 6 - length_digit, add.units = TRUE, cex = 0.6, tick.len = 3)

  ## gene track
  if (gene_track) {
    if (missing(TxDb)) {
      stop("TxDb is required. Please install the corresponding TxDb package")
    }
    genes.data <- karyoploteR::makeGenesDataFromTxDb(TxDb, karyoplot = kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data = genes.data, data.panel = 2, r0 = 0.7, r1 = 0.85, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
    karyoploteR::kpAddLabels(kp, labels = "Gene track", data.panel = 2, r0 = 0.7, r1 = 0.85, cex = 0.6, col = "black")
  }

  ## peak track
  if (peak_track) {
    if (missing(peaks)) {
      stop("peaks are required")
    }
    peaks_range <- regioneR::toGRanges(peaks)
    # karyoploteR::kpPlotRegions(kp, data = peaks_range, data.panel = 2, r0 = 0.9, r1 = 1, col = peaks_range$itemRgb)
    karyoploteR::kpPlotDensity(kp, data = peaks_range, data.panel = 2, r0 = 0.4, r1 = 0.1, window.size = 10^(length_digit - 1.1))
    karyoploteR::kpAddLabels(kp, labels = "Chromatin \n activity", data.panel = 2, r0 = 0.4, r1 = 0.1, cex = 0.6, col = "black")
  }

  ## plot peak
  ends <- regioneR::toGRanges(data_toplot[, c("Peak_chr", "Peak_start", "Peak_end")])
  karyoploteR::kpPlotRegions(kp, data = ends, r0 = 0.2, r1 = 0.3, col = "red")
  karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, x1 = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y0 = 0.5, y1 = 0.4, length = 0.05)
  karyoploteR::kpText(kp, labels = paste0("peak ", chr, ":", data_toplot$Peak_start, "-", data_toplot$Peak_end), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y = 0.65, col = "black", cex = 0.6)

  ### top TF(s)
  TF_binding_prob <- data_toplot[, "TFs.binding.prob."]
  TFs <- gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
  TFs <- TFs[TFs != ""]
  karyoploteR::kpText(kp, labels = paste0("top binding TF(s): ", paste(head(TFs, 3), collapse = " ")), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end) / 2, y = 0.55, col = "black", cex = 0.6)


  ## plot genes
  for (i in 1:nrow(data_toplot)) {
    starts <- regioneR::toGRanges(data_toplot[i, c("Gene_chr", "Gene_TSS", "Gene_TSS")])
    karyoploteR::kpPlotRegions(kp, data = starts, r0 = 0.2, r1 = 0.3, col = "blue")
    karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.3, r1 = 0.5, col = "grey")
    karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[, "Gene_TSS"], y = 0.05, col = "black", cex = 0.6)
    karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.1, y1 = 0.18, length = 0.05)
  }

  ## title
  bb <- karyoploteR::getMainTitleBoundingBox(kp)
  x <- (bb$x0 + bb$x1) / 2
  y <- 0.9

  graphics::text(x = x, y = y, labels = paste0("Circuits containing peak ", peak_chr, ":", peak_start, "-", peak_end), cex = 1, font = 2)
}
