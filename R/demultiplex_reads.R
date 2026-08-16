#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
demultiplex_reads <- function (reads, out.file = "None") 
{
  if (is.character(reads)) {
    info <- read.csv(out.file)
  }else {
    info <- reads
  }
  message("\n\n\n================================#| STEP 5: DEMULTIPLEXING |#================================")
  getmode <- function(v) {
    uniqv <- BiocGenerics::unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  #TODO these should both be inputs
  SampleIndex_WhiteList <- scMOSAIC::SampleIndex_WhiteList
  Multiome_key <- scMOSAIC::Multiome_key
  
  message(paste("number of reads: ", nrow(info)))
  info <- merge(info, SampleIndex_WhiteList[, c(1, 3)], by.x = "sampleid", 
                by.y = "index")
  info$barcode_sample <- paste(info$barcode, info$sample, sep = "_")
  message(paste("unique cells recovered: ", length(BiocGenerics::unique(info$barcode_sample))))
  
  message("\n\n\n================================#| STEP 5.1: REPEAT COUNTING |#================================")
  info$cag_repeats <- count_cag_repeats(info)
  
  cell_info <- parallel::mclapply(BiocGenerics::unique(info$barcode_sample), 
                                  FUN = function(cell) {
                                    sub <- info[info$barcode_sample == cell, ]
                                    res <- data.frame(atac = sub$barcode[1], 
                                                      cag_repeats_mean = mean(sub$cag_repeats), 
                                                      cag_repeats_median = median(sub$cag_repeats), 
                                                      cag_repeats_mode = getmode(sub$cag_repeats), 
                                                      cag_repeats_var = var(sub$cag_repeats), cag_repeats_sd = sd(sub$cag_repeats), 
                                                      number_of_reads = nrow(sub), sample = sub$sample[1], 
                                                      barcode_sample = cell)
                                    return(res)
                                  }, mc.cores = parallel::detectCores() - 1)
  cell_info <- do.call(rbind, cell_info)
  message("cell per sample recovered: ")
  message(BiocGenerics::table(cell_info$sample))
  cell_info <- merge(cell_info, Multiome_key, by = "atac")
  colSums(is.na(cell_info))
  message(paste("Number of unique cells: ", nrow(cell_info)))
  message(paste("Min number of reads per cell: ", min(cell_info$number_of_reads)))
  message(paste("Max number of reads per cell: ", max(cell_info$number_of_reads)))
  message(paste("Mean number of reads per cell: ", mean(cell_info$number_of_reads)))
  if (out.file == "None") {
    return(cell_info)
  }
  else {
    write.csv(cell_info, file = out.file)
    message("cell level info saved")
  }
}