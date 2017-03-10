
#' Parse Capture Hi-C data as \code{GInteractions} object.
#'
#' In the paper Mifsud et al. 2015
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/25938943}) promoter caputure Hi-C
#' data is provied for two cell types. This data can be downloaded from the
#' supplemental material of the studay as .zip archive:
#' \url{http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip}.
#' Each file in this zip archive can be parsed separately with this function as
#' an \code{\link[InteractioSet]{GInteractions}} object.
#'
#' @param infile Input file with interactions.
#' @param seqInfo A seqinfo object holding the length of each chromosome.
#' @return A \code{\link[InteractioSet]{GInteractions}} object
#' @export
parseCaptureHiC <- function(infile, seqInfo=NULL){

  # get classes of each column for faster parsing
  classes <- sapply(read.delim(infile,
                                nrows = 5,
                                header = TRUE,
                                stringsAsFactors = FALSE ), class)

  # parse input file as data frame
  df <- read.delim(infile,
                       header = TRUE,
                       colClasses = classes,
                       stringsAsFactors = FALSE)

  gr1 <- GenomicRanges::GRanges(
    df[,1],
      IRanges::IRanges(df[,2], df[,3]),
      seqinfo = seqInfo,
      symbol = df$Symbol,
      ENSG = df$Ensembl.Gene.ID,
      expresssion.quartile = df$expresssion.quartile)
  gr2 <- GenomicRanges::GRanges(
    df[,7],
    IRanges::IRanges(df[,8], df[,9]),
    seqinfo = seqInfo,
    symbol = df$Symbol.1,
    ENSG = df$Ensembl.Gene.ID.1,
    expresssion.quartile = df$expresssion.quartile.1)

  # construct GInteraction
  gi <- InteractionSet::GInteractions(
    gr1,
    gr2,
    raw.count = df$raw.count,
    log.observed.expected = df$log.observed.expected.
  )

  return(gi)
}
