#' Calculate Allele length 
#'
#' This function was designed to estimate the length of both the 'normal' inherited allele and the inherited expanded allele specifically in
#' CAG repeat disorders such as Huntingtons disease - but can be adapted for use in other expansion disorders.
#' This fuction expects a dataframe with 
#' 
#' 
#'
#' @param repeat_counts a numeric or integer vector corresponding to repeats 
#' @param scale if the data should be scaled; reccommended do not change
#' @param center if the data should be centered; reccommended do not change
#' @param plot provide some helpful plots to visualize what function is doing
#' @param z.norm.thresh the zscore threshold used to remove outliers in data. NOTE: long tails at end of expanded allele will skew expanded allele calling; recommended to remove outliers
#' @return A dataframe containing the insert, barcode ID, sample ID, and direction of the insert 
#' @export
calculate_alleles <- function(repeat_counts, scale = TRUE, center = TRUE, plot = FALSE, z.norm.thresh = 3){
  
  #take all the reads and scale them using a Z normalization; combine with original vector
  data <- data.frame(raw = repeat_counts, #input vector
                     scaled = scale(repeat_counts, #z normalize
                                           center = center, 
                                           scale = scale))
  
  #provide a plot to visualize the z.norm.thresh
  if(plot){
    
  p1 <- ggplot(data, aes(x = raw)) + 
    geom_histogram() + 
    ggtitle("raw data")
  
  p2 <- ggplot(data, aes(x = raw)) + 
    geom_histogram() + 
    ggtitle("z-normalized data")
  p3 <- ggplot(data[data$scaled < z.norm.thresh], aes(x = raw)) + 
    geom_histogram() + 
    ggtitle(paste0("z-normalized data within threshold (|value| < ", z.norm.thresh, ")"))
  
  print(p1 + p2 + p3)
  
  }
    
  #calculate modes assuming 2 modes; using only the data that falls withing 3 SD of the mean; ie remove outliers and calculate modes
  res <- locmodes(data = data[abs(repeat_counts_scaled) < 3,]$cag_repeats, #remove outliers
                  mod0 = 2, #assume two modes one for each allele
                  display = plot)
  
  #extract the values for each of mode and antimode
  locations <- round(res$locations)
  mode_info <- data.frame(wt = locations[1], anitmode = locations[2], hd = locations[3])
  
  return(mode_info)
}
