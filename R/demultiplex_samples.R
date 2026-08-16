#' Parse PacBio Reads 
#'
#' This function loads a fastq and parse the reads into its relevant poart. It assumes that the fastq is alread oriented from 5' prime to 3' prime
#' it searches for any gene between the exon1p1 and exon1p2 section
#' it extracts the insert, sample ID and barcode ID
#' 
#'
#' @param input Path to the input file
#' @param out.folder the location of the output file, where you want your results to be saved #TODO do we want this or should we just return the table
#' @param run.name *optional* for naming the output files
#' @param Read2N element from 10X
#' @param SampleIndexes optional used to subset and filter the reads based on if you know what the sample indexes you are expecting
#' @param Read2N.max.mm the max allowed mismatch for the read2n element
#' @param sampleid.length how long the 10 sample barcode is expected to be
#' @param min.reads the minimum number of reads needed from a sampleid to be considered non - junk
#' @param sampleid.max.mm the max mismatch allowed for barcode from the SampleIndexes to be considered a valid sampleID barcode
#' @param file.format if the file is in fastq or not, fasta?
#' @export
demultiplex_samples <- function(input, #List("data/lima.fl.5p--3p.fastq", "data/m54328U_241216_223053.lima.fl.5p--3p.fastq")
                        out.folder = "results/demultiplex_samples_testing/",#"results/all_data_htt_hits_info3.csv"
                        run.name = "test",
                        Read2N = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                        SampleIndexes = NULL, 
                        Read2N.max.mm = 1,
                        sampleid.length = 8,
                        sampleid.max.mm = 1,
                        min.reads = 1000,
                        file.format = "fastq", 
                        n.cores = parallel::detectCores()-1){
  
  #checks to make sure in and out files / paths exist
  if(!dir.exists(out.folder)){
    stop("out.folder does not exist.")
  }
  if(!stringr::str_detect(out.folder, "\\/$")){
    message("out.folder does not contain terminal / adding one")
    out.folder <- paste0(out.folder, "/")
  }
  if(!file.exists(input)){
    stop("input file does not exist.")
  }
  
  message("foo1")
  
  ###STEP 0 setup
  message("================================#| STEP 0: SETUP |#================================")
  
  #reduce valid sample index whitelist to only those we are expecting to find within our dataset
  #TODO make this optional in sampleIndex is null
  #TODO might want to make a different copy that contains all whitelist items in case user does a mistake and doesnt include an index 
  if(!is.null(SampleIndexes)){
    SampleIndex_WhiteList.filtered <- SampleIndex_WhiteList[SampleIndex_WhiteList$sample %in% SampleIndexes,]
    SampleIndex_WhiteList.error <- SampleIndex_WhiteList[!SampleIndex_WhiteList$sample %in% SampleIndexes,]
  }else{
    SampleIndex_WhiteList.filtered <- SampleIndex_WhiteList
    SampleIndex_WhiteList.error <- NULL
  }
  
  
  
  #set the working directory
  #setwd("/net/bmc-lab4/data/kellis/users/sfass/projects/HD")
  message(paste("working directory: ", getwd()))
  
  # -------------------------
  # Input handling
  # -------------------------
  if (inherits(input, "DNAStringSet")) {
    d <- input
    message("input is DNAStringSet")
  } else if (is.character(input) && length(input) == 1L) {
    message(paste0("reading: ", input))
    d <- Biostrings::readDNAStringSet(input, format = file.format)
    message("1 file loaded")
  } else if (is.character(input) && length(input) > 1L) {
    d <- do.call(
      c,
      lapply(input, function(file) {
        message(paste0("reading: ", file))
        Biostrings::readDNAStringSet(file, format = file.format)
      })
    )
    message(paste0(length(input), " files loaded"))
  } else {
    stop("input must be a DNAStringSet, a file path, or a character vector of file paths.")
  }
  
  if (length(d) == 0L) {
    stop("No reads loaded.")
  }
  
  if (is.null(names(d))) {
    names(d) <- paste0("read_", seq_along(d))
  }
  
  message("\n\n\n================================#| STEP 1: FILTERING |#================================")
  
  message(paste("Reads pre-filtering: ", length(d)))
  
  # Function to identify reads containing a specific sequence with mismatches
  #TODO how should I handle reads with multiple hits for an element?
  filter_reads <- function(target_seq, reads, max_mismatch) {
    hits <- Biostrings::vmatchPattern(target_seq, reads, max.mismatch = max_mismatch)
    #message(BiocGenerics::table(S4Vectors::elementNROWS(hits)))
    return(reads[S4Vectors::elementNROWS(hits) == 1])  # Subset reads containing the pattern
  }
  
  #filter the barcodes so each barcode can match a element in the whitelist with a max mismatch of one, 
  #only barcodes with a single match are retained and the whitelist element is returned
  #TODO for some reason this function is super memory heavy, probably because we are making an instance of whitelist for every barcode
  filter_barcodes <- function(barcodes, whitelist, max_distance = sampleid.max.mm) {
    #convert our whitelist to a dnastring set
    whitelist_set <- Biostrings::DNAStringSet(whitelist)
    
    #interate through each barcode and see if we have a match in the whitelist, within 1 mismatch range with only one unique match
    valid_barcodes <- parallel::mclapply(barcodes, FUN = function(barcode) {
      
      #look for barcodes that have a match in the whitelist with a max.mismatch of 1
      matches <- Biostrings::vcountPattern(pattern = as.character(barcode), whitelist_set, max.mismatch = max_distance)
      
      #if there are multiple matches or no matches to a given whitelist barcode then return NA, 
      #else return the corresponding barcode from the whitelist
      if (sum(matches) == 1) {
        return(whitelist[which(matches == 1)])
      }else{
        return(NA)
      }
    }, mc.cores = parallel::detectCores() %/% 2) #If we are going over out memory budget we can reduce this futher, 
    
    #BiocGenerics::unlist, remove barcodes without a close match, and convert back to desired format
    valid_barcodes <- BiocGenerics::unlist(valid_barcodes) #I guess the name sticks with the barcode
    valid_barcodes <- valid_barcodes[!is.na(valid_barcodes)]
    valid_barcodes <- Biostrings::DNAStringSet(valid_barcodes)
    
    return(valid_barcodes)
  }
  
  #TODO we can limit the search space to just the first 100bp
  
  #filter on r2 elements
  d <- filter_reads(Read2N, d, Read2N.max.mm)
  message(paste("Reads Read2N filtered: ", length(d)))
  
 
  ###STEP 2.4 - extract the sampleid
  message("\n\n\n================================#| STEP 2: SAMPLEID EXTRACTION |#================================")
  
  #remove reads where the end index of Read2N is equal to the width of the read, ie if any read ends with the read2n (ie no sampleid) then remove
  d <- d[!(BiocGenerics::unlist(Biostrings::endIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm))) + 1) > BiocGenerics::width(d)]
  
  #get the sample ID
  sampleids <- Biostrings::subseq(d, 
                                  start = BiocGenerics::unlist(Biostrings::endIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm))) + 1,
                                  end = BiocGenerics::width(d))
  message(paste0("number of unique sampleIDs found: ", length(unique(sampleids))))
  #message("2")
  
  #TODO can we leave the ones that are slightly off the expected length for recover in filter barcodes section?
  #make sure all the sampleids are the right length
  print(table(BiocGenerics::width(sampleids)))
  sampleids <- sampleids[BiocGenerics::width(sampleids) == sampleid.length]
  message(paste("correct length sampleids found: ", length(sampleids)))
  
  #message("3")
  
  #check for exact matches within the SampleIndex whitelist
  samp_exact_matches <- as.character(sampleids) %in% SampleIndex_WhiteList.filtered$index
  message(paste("sampleids recovered from SampleIndexes provided by user with exact matches: ", sum(samp_exact_matches)))
  
  if(!is.null(SampleIndex_WhiteList.error)){
  #lets make sure the user didnt F up and put in the wrong whitelists
  samp_error_matches <- as.character(sampleids) %in% SampleIndex_WhiteList.error$index
  message(paste("sampleids recovered from whitelist with exact matches: ", sum(samp_error_matches)))
  
  #if we observe a large amount of reads attributed to a sample not provided by user issue an alert
  sampleid.error.rate <- sum(samp_error_matches) / sum(samp_exact_matches)
  if(sampleid.error.rate > .05){
    message("WARNING possible mistake in user SampleIndexes")
    samp_error_matches_df <- as.data.frame(sampleids[samp_error_matches])
    samp_error_matches_df <- merge(samp_error_matches_df, SampleIndex_WhiteList.error, by.x = "x", by.y = "index", all.x = TRUE)
    message(paste(sampleid.error.rate, "% of sampleid founds matches to whitelist but not provided by user: did you forget to provide some SampleIndexes?"))
    et <- table(samp_error_matches_sampleids$sample)
    message(paste(names(et), ": ", et, "\n"))
  }
  }
  
  #message("4")
  
  #check for close sampleid index matches 
  samp_close_match <- filter_barcodes(sampleids[!samp_exact_matches], SampleIndex_WhiteList.filtered$index)
  message(paste("sampleids recovered from whitelist with close matches (1 mismatch, 1 unique barcode): ", length(samp_close_match)))
  
  #combine the close and exact matches
  sampleids <- c(samp_close_match, sampleids[samp_exact_matches])
  
  #info about the barcode filtering 
  message(paste("sampleids recovered from whitelist total: ", length(sampleids)))
  
  #message("5")
  
  #get unique SampleIndexes
  usid <- BiocGenerics::unique(sampleids)
  
  #make a unique SampleIndexes frequency table
  usidf <- BiocGenerics::table(sampleids)
  
  #output some information about sampleids
  message(paste("Unique SampleIndexes found: ", length(usid)))
  message(paste("SampleIndexes mean frequency:", mean(usidf)))
  message(paste("SampleIndexes median frequency:", median(usidf)))
  message(paste("SampleIndexes max frequency:", max(usidf)))
  
  #filter d down to only those with sampleids passing filtering
  d <- d[names(d) %in% names(sampleids)]
  
  message("\n\n\n================================#| STEP 2: fastq output |#================================")
  
  
  #lets make a table to decide on which indeces to save
  sampleids_df <- as.data.frame(sampleids)
  sampleids_df <- merge(sampleids_df, SampleIndex_WhiteList.filtered, by.x = "x", by.y = "index", all.x = TRUE)
 
  sidt <- table(sampleids_df$sample)
  message(paste(names(sidt), ": ", sidt, "\n"))
  
  SampleIndex_WhiteList.filtered <- SampleIndex_WhiteList.filtered[SampleIndex_WhiteList.filtered$sample %in% names(sidt[sidt > min.reads]),]
 
  #sample <- "A1"
  parallel::mclapply(unique(SampleIndex_WhiteList.filtered$sample), function(sample){
    
    #for each sampleid from the whitelist, gather all the indexes that match this ID
    sample_indeces <- SampleIndex_WhiteList[SampleIndex_WhiteList$sample == sample,]$index
    
    #subset all the reads to only those that have relevant sampleid indeces 
    d1 <- d[names(d) %in% names(sampleids[as.character(sampleids) %in% sample_indeces])]
    
    #write these to their own fastq file
    Biostrings::writeXStringSet(d1, filepath = paste0(out.folder, run.name, "_", sample, "_reads.fastq.gz"), format = "fastq", compress = TRUE)
  }, mc.cores = n.cores)
  
  message("wrote info")
  
}

