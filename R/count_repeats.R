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
    stopifnot(is.data.frame(i1))
    stopifnot(all(c("insert", "strand") %in% colnames(i1)))
    
    inserts <- as.character(i1$insert)
    strands <- as.character(i1$strand)
    
    count_one <- function(seq, strand, min.repeats) {
      if (is.na(seq) || is.na(strand) || !nzchar(seq)) {
        return(0L)
      }
      
      pattern <- if (strand == "fwd") "(?:CAG)+" else "(?:CTG)+"
      
      m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
      
      if (length(m) == 1L && m[1] == -1L) {
        return(0L)
      }
      
      match_lengths <- attr(m, "match.length")
      repeat_counts <- match_lengths / 3L
      repeat_counts <- repeat_counts[repeat_counts > min.repeats]
      
      if (length(repeat_counts) == 0L) {
        return(0L)
      }
      
      sum(repeat_counts)
    }
    
    out <- mapply(
      FUN = count_one,
      seq = inserts,
      strand = strands,
      MoreArgs = list(min.repeats = min.repeats),
      USE.NAMES = FALSE
    )
    
    names(out) <- rownames(i1)
    out
}






