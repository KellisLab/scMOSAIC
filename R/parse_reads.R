#' Parse PacBio Reads 
#'
#' This function loads a fastq and parse the reads into its relevant poart. It assumes that the fastq is alread oriented from 5' prime to 3' prime
#' it searches for any gene between the exon1p1 and exon1p2 section
#' it extracts the insert, sample ID and barcode ID
#' 
#'
#' @param input Path to the input file
#' @param out.file the location of the output file, where you want your results to be saved #TODO do we want this or should we just return the table
#' @param plot boolean of if you want plots to be created along the way to visualize the process
#' @param P5element TODO is this still necessary?
#' @param Spacer spacer element from 10X
#' @param Read1N element from 10X
#' @param Read2N element from 10X
#' @param P7Element element from 10X
#' @param HTTexon1_p1 TODO change this name to something more general; desired sequence match start sequence
#' @param HTTexon1_p2 TODO change this name to something more general; desired sequence match start sequence
#' @param SampleIndexes optional used to subset and filter the reads based on if you know what the sample indexes you are expecting
#' @param P5element.max.mm the max allowed mismatch for the P5 element TODO still necessary?
#' @param Spacer.max.mm the max allowed mismatch for the P5 element TODO still necessary?
#' @param Read1N.max.mm the max allowed mismatch for the read1n element
#' @param Read2N.max.mm the max allowed mismatch for the read2n element
#' @param P7Element.max.mm the max allowed mismatch for the TODO still necessary?
#' @param HTTexon1_p1.max.mm the max allowed mismatch for the first part of the desired sequence matched gene #TODO make this not names HTT
#' @param barcode.length how long the 10X barcodes are expected to be
#' @param barcode.max.mm the max mismatch allowed for barcode from the whitelist to be considered a valid barcode
#' @param sampleid.length how long the 10 sample barcode is expected to be
#' @param sampleid.max.mm the max mismatch allowed for barcode from the SampleIndexes to be considered a valid sampleID barcode
#' @param file.format if the file is in fastq or not, fasta?
#' @return A dataframe containing the insert, barcode ID, sample ID, and direction of the insert 
#' @export
parse_reads <- function(input, #List("data/lima.fl.5p--3p.fastq", "data/m54328U_241216_223053.lima.fl.5p--3p.fastq")
                        out.file,#"results/all_data_htt_hits_info3.csv"
                        plot = FALSE, 
                        P5element = "AATGATACGGCGACCACCGAGATCTACAC",
                        Spacer = "CGCGTCTG",
                        Read1N = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                        Read2N = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                        P7Element = "ATCTCGTATGCCGTCTTCTGCTTG",
                        HTTexon1_p1 = "ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC",
                        HTTexon1_p2 = "CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCGGCCCGGCTGTGGCTGAGGAGCCGCTGCACCGACC",
                        SampleIndexes = c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), 
                        P5element.max.mm = 1,
                        Spacer.max.mm = 0,
                        Read1N.max.mm = 1,
                        Read2N.max.mm = 1,
                        P7Element.max.mm = 1,
                        HTTexon1_p1.max.mm = 3,
                        HTTexon1_p2.max.mm = 6, 
                        barcode.length = 16,
                        barcode.max.mm = 1,
                        sampleid.length = 8,
                        sampleid.max.mm = 1,
                        file.format = "fastq"){
  
  #TODO put in checks to make sure in and out files / paths exist
  message("foo")
  
  ###STEP 0 setup
  message("================================#| STEP 0: SETUP |#================================")
  
  #declare some sequences that are important structures in the reads
  HTTexon1_p1_rvcmpt <- as.character( Biostrings::reverseComplement(Biostrings::DNAString(HTTexon1_p1)))
  HTTexon1_p2_rvcmpt <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(HTTexon1_p2)))
  
  #reduce valid sample index whitelist to only those we are expecting to find within our dataset
  #TODO make this optional in sampleIndex is null
  #TODO this should be an input
  SampleIndex_WhiteList <- scMOSAIC::SampleIndex_WhiteList
  SampleIndex_WhiteList <- SampleIndex_WhiteList[SampleIndex_WhiteList$sample %in% SampleIndexes,]
  
  #TODO this should be an input
  Barcode10X_WhiteList <- scMOSAIC::Barcode10X_WhiteList
  
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
  
  #TODO we can limit the search space to just the first 100bp
  #filter on spacer elements
  d <- filter_reads(Spacer, d, Spacer.max.mm)
  message(paste("Reads Spacer filtered: ", length(d)))
  
  #filter on r1 elements
  d <- filter_reads(Read1N, d, Read1N.max.mm)
  message(paste("Reads Read1N filtered: ", length(d)))
  
  #filter on r2 elements
  d <- filter_reads(Read2N, d, Read2N.max.mm)
  message(paste("Reads Read2N filtered: ", length(d)))
  
  #filter the reads down so only reads where the read2n is 8 away from the end of the sample remain
  d <- d[BiocGenerics::width(d) - (BiocGenerics::unlist(Biostrings::endIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm)))) == sampleid.length]
  message(paste("Reads SampleIndexN position filtered: ", length(d)))
  
  #TODO we can limit the search space to just the last 100bp?
  #remove all reads that dont have a spacer starting at location 17
  #TODO should we do spacer plus read1?
  d <-  d[BiocGenerics::unlist(Biostrings::startIndex(Biostrings::vmatchPattern(Spacer, d, max.mismatch = Spacer.max.mm))) == barcode.length + 1]
  message(paste("Reads 10XBarcode position filtered: ", length(d)))
  
  #check for HTT gene accounting for both forward and reverse complement strands
  #(p1 & p2)
  p1_matches <- S4Vectors::elementNROWS(Biostrings::vmatchPattern(HTTexon1_p1, d, max.mismatch = HTTexon1_p1.max.mm)) == 1 #p1 matches
  message(paste("number of forward HTT p1 containing reads found: ", sum(p1_matches)))
  
  p2_matches <- S4Vectors::elementNROWS(Biostrings::vmatchPattern(HTTexon1_p2, d, max.mismatch = HTTexon1_p2.max.mm)) == 1 #p2 matches
  message(paste("number of forward HTT p2 containing reads found: ", sum(p2_matches)))
  
  fwd <- (p1_matches & p2_matches)
  message(paste("number of forward HTT (p1 & p2) reads found: ", sum(fwd)))
  
  #add metadata indicating if each read was forward
  S4Vectors::mcols(d)$fwd <- fwd
  
  #find reads that contain reverse complement p1 and p2
  p1_matches_rvcmpt <- S4Vectors::elementNROWS(Biostrings::vmatchPattern(HTTexon1_p1_rvcmpt, d, max.mismatch = HTTexon1_p1.max.mm)) == 1
  message(paste("number of reverse complement HTT p1 containing reads found: ", sum(p1_matches_rvcmpt)))
  
  p2_matches_rvcmpt <- S4Vectors::elementNROWS(Biostrings::vmatchPattern(HTTexon1_p2_rvcmpt, d, max.mismatch = HTTexon1_p2.max.mm)) == 1
  message(paste("number of reverse complement HTT p2 containing reads found: ", sum(p2_matches_rvcmpt)))
  
  rvcmpt <- (p1_matches_rvcmpt & p2_matches_rvcmpt)
  message(paste("number of reverse complement HTT reads found: ", sum(rvcmpt)))
  
  #add metadata indicating if each read was reverse complement
  S4Vectors::mcols(d)$rvcmpt <- rvcmpt
  
  #include forward and reverse complement HTT reads, toss the rest
  d <- d[fwd | rvcmpt]
  message(paste("Reads HTTexon1_p2 & HTTexon1_p1 filtered: ", length(d)))
  
  #make a better metadata column for later use
  #TODO make this more elegant
  S4Vectors::mcols(d)$strand <- c("foo") #make a dummy column to populate
  d[d@elementMetadata$fwd]@elementMetadata$strand <- c("fwd")
  d[d@elementMetadata$rvcmpt]@elementMetadata$strand <- c("rvcmpt")
  
  #TODO can we print this table off?
  #message("table of reads captured in forward and reverse complement confirmations:")
  #BiocGenerics::table(d@elementMetadata$strand)
  
  ###STEP 2 - extract the important information, aka barcode, insert and sampleID###
  #message("================================#| STEP 2: EXTRACTION |#================================")
  
  
  ###STEP 2.1 - extract barcodes, filter data by barcode length, filter data by removing reads with garbage barcodes 
  message("\n\n\n================================#| STEP 2.1: BARCODE EXTRACTION |#================================")
  
  #extract the barcodes, I guess these should always start at the first position and be exactly the length of the barcode as specificed by the user 
  barcodes <- Biostrings::subseq(d, start = 1, end = barcode.length)
  
  #2.3 filter the barcodes
  
  #filter the barcodes so each barcode can match a element in the whitelist with a max mismatch of one, 
  #only barcodes with a single match are retained and the whitelist element is returned
  #TODO for some reason this function is super memory heavy, probably because we are making an instance of whitelist for every barcode
  filter_barcodes <- function(barcodes, whitelist, max_distance = barcode.max.mm) {
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
  
  #Get all the exact matches between the barcodes and the whitelist
  exact_matches <- as.character(barcodes) %in% Barcode10X_WhiteList
  message(paste("Barcodes recovered from whitelist with exact matches: ", sum(exact_matches))) #add in to the dataframe .. etc
  
  #get all the matches that are not exact, that are within 1 mismatch, but only have 1 unique match in the whitelist
  close_matches <- filter_barcodes(barcodes[!exact_matches], Barcode10X_WhiteList)
  message(paste("Barcodes recovered from whitelist with close matches (", barcode.max.mm, " mismatch, 1 unique barcode): ", length(close_matches)))
  
  #combine the close matches and the exact matches
  barcodes <- c(close_matches, barcodes[exact_matches])
  
  #info about the barcode filtering 
  message(paste("Barcodes recovered from whitelist total: ", length(barcodes)))
  
  #get unique barcodes
  ubc <- BiocGenerics::unique(barcodes)
  
  #make a unique barcode frequency table
  bcf <- BiocGenerics::table(barcodes)
  
  #output some information about barcodes
  message(paste("Unique barcodes found: ", length(ubc)))
  message(paste("Barcodes mean frequency:", mean(bcf)))
  message(paste("Barcodes median frequency:", median(bcf)))
  message(paste("Barcodes max frequency:", max(bcf)))
  
  #2.2 optionally visualize the similarity of the barcodes
  
  #Plotting barcode similarity#
  if(plot){
    
    #TODO improve and add titles    
    graphics::hist(bcf)
    graphics::hist(log(bcf))
    
    if(length(ubs) > 10000){
      message("too many unique barcodes to visualize, subsetting...")
      ubc <- sample(ubc, size = 10000, replace = FALSE)
    }
    
    #
    dist_matrix <- as.matrix(stringDist(ubc, method = "hamming"))
    
    # Function to group barcodes by distance threshold
    get_unique_count <- function(dist_matrix, threshold) {
      groups <- stats::hclust(stats::as.dist(dist_matrix), method = "single") %>%
        stats::cutree(h = threshold)
      length(BiocGenerics::unique(groups))
    }
    
    # Compute unique counts for different distance thresholds
    thresholds <- seq(0, max(dist_matrix), by = 1)
    unique_counts <- BiocGenerics::sapply(thresholds, function(th) get_unique_count(dist_matrix, th))
    
    # Prepare data for plotting
    data <- data.frame(
      Threshold = thresholds,
      UniqueBarcodes = unique_counts
    )
    
    # Plot the results
    ggplot2::ggplot(data, aes(x = Threshold, y = UniqueBarcodes)) +
      geom_line(size = 1, color = "blue") +
      geom_point(size = 3, color = "red") +
      labs(
        title = "Number of Unique Barcodes by Distance Threshold",
        x = "Distance Threshold (Hamming Distance)",
        y = "Unique Barcodes"
      ) +
      theme_minimal()
  }
  
  #filter d down to only those with barcodes passing filtering
  d <- d[names(d) %in% names(barcodes)]
  message(paste("barcode whitelist filtered reads remaining: ", length(d)))
  
  
  ###STEP 2.4 - extract the sampleid
  message("\n\n\n================================#| STEP 2.3: SAMPLEID EXTRACTION |#================================")
  
  #message("1")
  #get the sample ID
  sampleids <- Biostrings::subseq(d, 
                                  start = BiocGenerics::unlist(Biostrings::endIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm))) + 1,
                                  end = BiocGenerics::width(d))
  message(paste0("number of unique sampleIDs found: ", length(unique(sampleids))))
  #message("2")
  
  #make sure all the sampleids are the right length
  stopifnot(all(BiocGenerics::width(sampleids) == sampleid.length)) 
  
  #message("3")
  
  #check for exact matches within the SampleIndex whitelist
  samp_exact_matches <- as.character(sampleids) %in% SampleIndex_WhiteList$index
  message(paste("sampleids recovered from whitelist with exact matches: ", sum(samp_exact_matches)))
  
  #message("4")
  
  #check for close sampleid index matches 
  samp_close_match <- filter_barcodes(sampleids[!samp_exact_matches], SampleIndex_WhiteList$index)
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
  
  #filter d down to only those with barcodes passing filtering
  d <- d[names(d) %in% names(sampleids)]
  
  #filter down the cell barcodes to only those that also have a valid sampleid...
  barcodes <- barcodes[names(barcodes) %in% names(sampleids)]
  message(paste("SampleIndexes Filtered: ", length(d)))
  
  if(plot){
    
    #sampleid frequency histograms  
    graphics::hist(usidf)
    graphics::hist(log(usidf))
    
    dist_matrix <- as.matrix(stringDist(BiocGenerics::unique(sampleids), method = "hamming"))
    #dist_matrix <- as.matrix(stringDist(test, method = "levenshtein"))
    
    # Function to group barcodes by distance threshold
    get_unique_count <- function(dist_matrix, threshold) {
      groups <- stats::hclust(stats::as.dist(dist_matrix), method = "single") %>%
        stats::cutree(h = threshold)
      length(BiocGenerics::unique(groups))
    }
    
    # Compute unique counts for different distance thresholds
    thresholds <- seq(0, max(dist_matrix), by = 1)
    unique_counts <- BiocGenerics::sapply(thresholds, function(th) get_unique_count(dist_matrix, th))
    
    # Prepare data for plotting
    data <- data.frame(
      Threshold = thresholds,
      UniqueBarcodes = unique_counts
    )
    
    # Plot the results
    ggplot2::ggplot(data, aes(x = Threshold, y = UniqueBarcodes)) +
      geom_line(size = 1, color = "blue") +
      geom_point(size = 3, color = "red") +
      labs(
        title = "Number of Unique sampleIDs by Distance Threshold",
        x = "Distance Threshold",
        y = "Unique Barcodes"
      ) +
      theme_minimal()
  }
  
  ###STEP 2.3 - extract the inserts, TODO do I need to do any filtering here? should I do alignment here? 
  message("================================#| STEP 2.3: INSERT EXTRACTION |#================================")
  
  #extract the inserts 
  #TODO this if wrong! This should be the end index of spacer+read1n and start index of read2n
  # inserts <- Biostrings::subseq(d, 
  #                               start = BiocGenerics::unlist(startIndex(Biostrings::vmatchPattern(Read1N, d, max.mismatch = Read1N.max.mm))),
  #                               end = BiocGenerics::unlist(Biostrings::endIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm))))
  
  #new code
  inserts <- Biostrings::subseq(d, 
                    start = unlist(Biostrings::endIndex(Biostrings::vmatchPattern(paste0(Spacer,Read1N), d, max.mismatch = Spacer.max.mm + Read1N.max.mm))) + 1,
                    end = unlist(Biostrings::startIndex(Biostrings::vmatchPattern(Read2N, d, max.mismatch = Read2N.max.mm))) - 1)
  
  message(paste("Insert mean BiocGenerics::width:", mean(BiocGenerics::width(inserts))))
  message(paste("Insert median BiocGenerics::width:", median(BiocGenerics::width(inserts))))
  message(paste("Insert min BiocGenerics::width:", min(BiocGenerics::width(inserts))))
  message(paste("Insert max BiocGenerics::width:", max(BiocGenerics::width(inserts))))
  
  #inserts <- inserts[names(inserts) %in% names(sampleids)]
  
  
  ###STEP 3 - "demultiplex" the samples 
  ###STEP 3 - create CSV with barcode, insert, sample ID, gene name, insert length, CAG repeat value
  #TODO this isnt the demultiplexing step
  message("\n\n\n================================#| STEP 3: WRANGLE |#================================")
  
  # Convert to dataframes explicitly (ensure rownames exist)
  barcodes.df <- data.frame(barcode = as.character(barcodes))
  inserts.df <- data.frame(insert = as.character(inserts))
  sampleids.df <- data.frame(sampleid = as.character(sampleids))
  strand.df <- data.frame(strand = as.character(d@elementMetadata$strand))
  
  #ensure there are no duplicate names
  names(barcodes) <- make.unique(names(barcodes), sep = ".")
  names(inserts) <- make.unique(names(inserts), sep = ".")
  names(sampleids) <- make.unique(names(sampleids), sep = ".")
  names(d) <- make.unique(names(d), sep = ".")
  
  # Assign rownames explicitly (e.g., from original names of sequences)
  rownames(barcodes.df) <- names(barcodes)
  rownames(inserts.df) <- names(inserts)
  rownames(sampleids.df) <- names(sampleids)
  rownames(strand.df) <- names(d)
  
  # Reorder based on rownames while preserving rownames
  barcodes.df <- barcodes.df[order(rownames(barcodes.df)), , drop = FALSE]
  inserts.df <- inserts.df[order(rownames(inserts.df)), , drop = FALSE]
  sampleids.df <- sampleids.df[order(rownames(sampleids.df)), , drop = FALSE]
  strand.df <- strand.df[order(rownames(strand.df)), , drop = FALSE]
  
  #check to make sure all the names are lined up correctly 
  stopifnot(all(rownames(inserts.df) == rownames(barcodes.df)))
  stopifnot(all(rownames(inserts.df) == rownames(sampleids.df)))
  stopifnot(all(rownames(sampleids.df) == rownames(barcodes.df)))
  stopifnot(all(rownames(sampleids.df) == rownames(strand.df)))
  
  #combine all the information
  info <- cbind(barcodes.df, inserts.df, sampleids.df, strand.df)
  
  #write to outfile
  #TODO if print...
  if(out.file != "None"){
    utils::write.csv(info, file = out.file)
    message("wrote info to csv")
  }
  
  return(info)

}

