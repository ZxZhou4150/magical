#' Plot circuit(s)
#'
#' This function plots one circuit by its index.
#'
#' @seealso [plot_circuits_with_gene()] and [plot_circuits_with_peak()]
#'
#' @param data The data frame written by [MAGICAL_circuits_output()].
#' @param idx The index of the circuit to be plotted.
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `FALSE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' #library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #plot_circuits_with_idx(data, 1, T, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' @export

plot_circuits_with_idx = function(data, idx, gene_track=F, TxDb){
  if(length(idx)!=1){stop("idx should be a single number")}

  data_toplot = data[idx,]

  ## finding region
  chr = data_toplot$Gene_chr
  maxsite = max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite = min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit = floor(log10(maxsite - minsite))
  start = minsite - minsite %% 10^length_digit - 10^length_digit/2
  end = maxsite + (-maxsite %% 10^length_digit) + 10^length_digit/2

  region = regioneR::toGRanges(paste0(chr, ":", start, "-", end))
  kp = karyoploteR::plotKaryotype(zoom = region)

  ## gene track
  if(gene_track){
    if(missing(TxDb)){stop("TxDb is required. Please install the corresponding TxDb package")}
    genes.data = karyoploteR::makeGenesDataFromTxDb(TxDb,karyoplot=kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data=genes.data, data.panel = 2, r0 = 0.65, r1 = 0.8, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
  }

  ## finding start and end

  starts = regioneR::toGRanges(data_toplot[,c("Gene_chr","Gene_TSS","Gene_TSS")])
  ends = regioneR::toGRanges(data_toplot[,c("Peak_chr","Peak_start","Peak_end")])

  ## plotting
  karyoploteR::kpPlotRegions(kp, data = starts, r0=0, r1 = 0.1, col = "blue")
  karyoploteR::kpPlotRegions(kp, data = ends, r0=0, r1 = 0.1, col = "red")
  karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.1, r1 = 0.3, col = "grey")
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit-1), digit = 6-length_digit, add.units = TRUE, cex=0.6, tick.len = 3)
  karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[,"Gene_TSS"], y=0.5, col = "black", cex=0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.4, y1 = 0.2, length = 0.05)
  karyoploteR::kpText(kp, labels = paste0("peak ", chr,":",data_toplot$Peak_start, "-", data_toplot$Peak_end), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.55, col = "black", cex=0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot$Peak_start + data_toplot$Peak_end)/2, x1 = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y0 = 0.4, y1 = 0.2, length = 0.05)

  ### top TF(s)
  TF_binding_prob = data_toplot[,"TFs.binding.prob."]
  TFs = gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
  TFs = TFs[TFs!=""]
  if(length(TFs)>=3){
    karyoploteR::kpText(kp, labels = paste0("top TF(s): ", paste(head(TFs, 3), collapse = " ")) , chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.45, col = "black", cex=0.6)
  }
  else{karyoploteR::kpText(kp, labels = paste0("top TF(s): ",paste(TFs, collapse = " ")) , chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.45, col = "black", cex=0.6)}

  ### title
  bb = karyoploteR::getMainTitleBoundingBox(kp)
  x = (bb$x0+bb$x1)/2
  if(gene_track){y = 0.9}
  else{y = 0.8}

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
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `FALSE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' #library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #plot_circuits_with_gene(data, "ACOT11", T, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' @export

plot_circuits_with_gene = function(data, gene, gene_track=F, TxDb){
  if(length(gene)!=1){stop("gene should be a single gene symbol")}
  idx = which(data$Gene_symbol == gene)
  if(length(idx)==0){stop("gene not found")}

  data_toplot = data[idx,]

  ## finding region
  chr = data_toplot[1,"Gene_chr"]
  maxsite = max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite = min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit = floor(log10(maxsite - minsite))
  start = minsite - minsite %% 10^length_digit - 10^length_digit/2
  end = maxsite + (-maxsite %% 10^length_digit) + 10^length_digit/2

  region = regioneR::toGRanges(paste0(chr, ":", start, "-", end))

  kp = karyoploteR::plotKaryotype(zoom = region)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit-1), digit = 6-length_digit,add.units = TRUE, cex=0.6, tick.len = 3)

  ## gene track
  if(gene_track){
    if(missing(TxDb)){stop("TxDb is required. Please install the corresponding TxDb package")}
    genes.data = karyoploteR::makeGenesDataFromTxDb(TxDb,karyoplot=kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data=genes.data, data.panel = 2, r0 = 0.65, r1 = 0.8, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
  }

  ## plot gene
  starts = regioneR::toGRanges(data_toplot[1,c("Gene_chr","Gene_TSS","Gene_TSS")])
  karyoploteR::kpPlotRegions(kp, data = starts, r0=0, r1 = 0.1, col = "blue")
  karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[,"Gene_TSS"], y=0.6, col = "black", cex=0.6)
  karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.4, y1 = 0.2, length = 0.05)

  ## plot peaks
  for(i in 1:nrow(data_toplot)){
    ends = regioneR::toGRanges(data_toplot[i,c("Peak_chr","Peak_start","Peak_end")])
    karyoploteR::kpPlotRegions(kp, data = ends, r0=0, r1 = 0.1, col = "red")
    karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.1, r1 = 0.3, col = "grey")
    karyoploteR::kpText(kp, labels = paste0("peak #",i, " ", chr, ":", data_toplot[i,"Peak_start"] ,"-", data_toplot[i,"Peak_end"]), chr = chr, x = (data_toplot[i,"Peak_start"] + data_toplot[i,"Peak_end"])/2, y = 0.55, col = "black", cex=0.6)
    karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot[i,"Peak_start"] + data_toplot[i,"Peak_end"])/2, x1 = (data_toplot[i,"Peak_start"] + data_toplot[i,"Peak_end"])/2, y0 = 0.4, y1 = 0.2, length = 0.05)

    ### top TF(s)
    TF_binding_prob = data_toplot[i,"TFs.binding.prob."]
    TFs = gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
    TFs = TFs[TFs!=""]
    if(length(TFs)>=3){
      karyoploteR::kpText(kp, labels = paste0("top TF(s): ", paste(head(TFs, 3), collapse = " ")) , chr = chr, x = (data_toplot[i,"Peak_start"] + data_toplot[i,"Peak_end"])/2, y = 0.45, col = "black", cex=0.6)
    }
    else{karyoploteR::kpText(kp, labels = paste0("top TF(s): ",paste(TFs, collapse = " ")) , chr = chr, x = (data_toplot[i,"Peak_start"] + data_toplot[i,"Peak_end"])/2, y = 0.45, col = "black", cex=0.6)}
  }

  ### title
  bb = karyoploteR::getMainTitleBoundingBox(kp)
  x = (bb$x0+bb$x1)/2
  if(gene_track){y = 0.9}
  else{y = 0.8}

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
#' @param gene_track Whether to plot gene tracks around the region or not. Default is `FALSE`.
#' @param TxDb The TxDb object corresponding to the reference genome. Required if `gene_track` is `TRUE`.
#'
#' @return A karyoploteR plot.
#'
#' @examples
#' #library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' #data = read.table("MAGICAL_selected_regulatory_circuits.txt", header = T, sep = "\t")
#' #plot_circuits_with_peak(data,"chr4", 184817550, 184818412, T, TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' @export

plot_circuits_with_peak = function(data, peak_chr, peak_start, peak_end, gene_track=F, TxDb){
  idx = which(data$Peak_chr == peak_chr & data$Peak_start == peak_start & data$Peak_end == peak_end)
  if(length(idx)==0){stop("peak not found")}

  data_toplot = data[idx,]

  ## finding region
  chr = peak_chr
  maxsite = max(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  minsite = min(data_toplot$Gene_TSS, data_toplot$Peak_start, data_toplot$Peak_end)
  length_digit = floor(log10(maxsite - minsite))
  start = minsite - minsite %% 10^length_digit - 10^length_digit/2
  end = maxsite + (-maxsite %% 10^length_digit) + 10^length_digit/2

  region = regioneR::toGRanges(paste0(chr, ":", start, "-", end))

  kp = karyoploteR::plotKaryotype(zoom = region)
  karyoploteR::kpAddBaseNumbers(kp, tick.dist = 10^length_digit, minor.tick.dist = 10^(length_digit-1), digit = 6-length_digit,add.units = TRUE, cex=0.6, tick.len = 3)

  ## gene track
  if(gene_track){
    if(missing(TxDb)){stop("TxDb is required. Please install the corresponding TxDb package")}
    genes.data = karyoploteR::makeGenesDataFromTxDb(TxDb,karyoplot=kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
    genes.data <- karyoploteR::addGeneNames(genes.data)
    genes.data <- karyoploteR::mergeTranscripts(genes.data)

    karyoploteR::kpPlotGenes(kp, data=genes.data, data.panel = 2, r0 = 0.65, r1 = 0.8, gene.name.cex = 0.6, gene.name.position = "top", transcript.name.position = "top", add.gene.names = T, col = "grey")
  }

  ## plot peak
  ends = regioneR::toGRanges(data_toplot[,c("Peak_chr","Peak_start","Peak_end")])
  karyoploteR::kpPlotRegions(kp, data = ends, r0=0, r1 = 0.1, col = "red")
  karyoploteR::kpArrows(kp, chr = chr, x0 = (data_toplot$Peak_start + data_toplot$Peak_end)/2, x1 = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y0 = 0.4, y1 = 0.2, length = 0.05)
  karyoploteR::kpText(kp, labels = paste0("peak ", chr,":",data_toplot$Peak_start, "-", data_toplot$Peak_end), chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.55, col = "black", cex=0.6)

  ### top TF(s)
  TF_binding_prob = data_toplot[,"TFs.binding.prob."]
  TFs = gsub("\\(.*?\\)", "", strsplit(TF_binding_prob, ",\\s*")[[1]])
  TFs = TFs[TFs!=""]
  if(length(TFs)>=3){
    karyoploteR::kpText(kp, labels = paste0("top TF(s): ", paste(head(TFs, 3), collapse = " ")) , chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.45, col = "black", cex=0.6)
  }
  else{karyoploteR::kpText(kp, labels = paste0("top TF(s): ",paste(TFs, collapse = " ")) , chr = chr, x = (data_toplot$Peak_start + data_toplot$Peak_end)/2, y = 0.45, col = "black", cex=0.6)}

  ## plot genes
  for(i in 1:nrow(data_toplot)){
    starts = regioneR::toGRanges(data_toplot[i,c("Gene_chr","Gene_TSS","Gene_TSS")])
    karyoploteR::kpPlotRegions(kp, data = starts, r0=0, r1 = 0.1, col = "blue")
    karyoploteR::kpPlotLinks(kp, data = starts, data2 = ends, r0 = 0.1, r1 = 0.3, col = "grey")
    karyoploteR::kpText(kp, labels = data_toplot$Gene_symbol, chr = chr, x = data_toplot[,"Gene_TSS"], y=0.6, col = "black", cex=0.6)
    karyoploteR::kpArrows(kp, chr = chr, x0 = data_toplot$Gene_TSS, x1 = data_toplot$Gene_TSS, y0 = 0.4, y1 = 0.2, length = 0.05)
  }

  ## title
  bb = karyoploteR::getMainTitleBoundingBox(kp)
  x = (bb$x0+bb$x1)/2
  if(gene_track){y = 0.9}
  else{y = 0.8}

  graphics::text(x = x, y = y, labels = paste0("Circuits containing peak ", peak_chr, ":", peak_start, "-", peak_end), cex = 1, font = 2)
}
