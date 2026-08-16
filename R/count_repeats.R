#' count repeats
#'
#' This function 
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#' TODO fix so that the pattern is dynamic 
#'
#' @param i1 a dataframe containing a column with the name "insert" which has a character array of DNA
#' @param min.repeats is the number of consecutive repeats of a sequence needed to declare it a repeat
#' @return the repeated counts for the given sequence
#' @export
count_cag_repeats <- function(i1, min.repeats = 3) {
  # Function to count CAG repeats in a single sequence
  count_single_sequence <- function(rowname, i = i1) {
    
    row = i[rownames(i) == rowname,]
    
    # Regex pattern for the longest uninterrupted block of CAG repeats
    if(row$strand == "fwd"){
      pattern <- paste0("(CAG)+")
    }else{
      pattern <- paste0("(CTG)+")
    }
    
    # Find all matches of the pattern
    matches <- regmatches(row$insert, gregexpr(pattern, row$insert, perl = TRUE))[[1]]
    
    # Calculate the number of repeats for each match
    repeat_counts <- BiocGenerics::sapply(matches, function(match) nchar(match) / 3)
    
    #assuming that greater 3 is part of the CAG repeat region
    repeat_counts <- repeat_counts[repeat_counts > min.repeats]
    
    # Return the maximum number of repeats
    return(sum(repeat_counts))
  }
  
  # Apply the function in parallel using parallel::mclapply
  result <- parallel::mclapply(rownames(i1), count_single_sequence, mc.cores = parallel::detectCores() - 1)
  
  # Simplify the result to a vector
  result <- BiocGenerics::unlist(result)
  
  # Return a named vector where names are the sequences and values are repeat counts
  return(result)
}





