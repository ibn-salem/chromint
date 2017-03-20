

#' Get bins of a single chromosome.
#'
#' This function returns a \code{\link[GenomicRanges]{GRanges}} object with bins of a single
#' chromosome by given bin size.
#' @param chr Chromosome name
#' @param resolution The bin size resolution in base pairs
#' @param seqInfo A seqinfo object holding the length of each chromosome.
#' @return A GRanges object with bins
#' @keywords internal
#'
getBinGR <- function(chr, resolution, seqInfo){

  chrLen = GenomeInfoDb::seqlengths(seqInfo)[chr]
  starts = seq(1, chrLen, resolution)
  ends = starts+resolution-1
  ends = ifelse(ends < chrLen, ends, chrLen)
  n = length(starts)
  gr = GenomicRanges::GRanges(rep(chr, n),
    IRanges::IRanges(starts, ends),
    names=paste0("bin_", starts-1),
    seqinfo=seqInfo
  )

  return(gr)
}


#' Get bin offset for each chromosome.
#'
#' This function creates offset of bin indices for each chromosome.
#' @param chromosomes Character vecotr with chromosome identifiers.
#' @param resolution The bin size resolution in base pairs
#' @param seqInfo A seqinfo object holding the length of each chromosome.
#' @return A numeric vector with offset for each chromosome.
#' @keywords internal
getChrToBinOffset <- function(chromosomes, resolution, seqInfo){

	# get length of all chromosomes
	chrLen = GenomeInfoDb::seqlengths(seqInfo)[chromosomes]

	# get number of bins for each chromosomes
	nBins <- ceiling(chrLen / resolution)

	# get cummulative sum for offset
	binOffset <- cumsum(c(0, nBins))[1:length(chromosomes)]

	# fix names because of offset
	names(binOffset) <- chromosomes

  return(binOffset)
}


#' Parse Hi-C interactions from Rao et. al 2014 as
#' \code{\link[InteractionSet]{GInteractions}} object.
#'
#' This function parses interactions from a single experiment by a given
#' resolution. It takes a path to a directory with uncompressed raw data and
#' subdirectories for each cell-type. The data can be manually downloaded and
#' from this url:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525}.
#'
#' @param cell The cell type of the input experiment
#' @param resolution The matrix resolution of the Hi-C map in bp
#' @param baseDir The path to the directory with the Hi-C data.
#' @param seqInfo A seqinfo object of the genome.
#' @param interChromosomal logical idicating whether inter-chromosomal contacts
#'   should be parsed.
#' @param normalizeByExpected Should contact frequencies be normalized by
#'   distance. Default is FALSE.
#' @param mapQual String representing the mapping quality threshold. Defualt
#'   is "MAPQGE30".
#' @param norm String representation of the normalization method. Default is
#'   "KR". This is not implemented yet.
#' @keywords Hi-C
#' @return \code{\link[InteractionSet]{GInteractions}} object with all cis- and
#'   trans-interactions
#' @export
parseRaoHiCtoGI <- function(cell, resolution, baseDir, seqInfo,
                            interChromosomal = TRUE,
                            normalizeByExpected = FALSE,
                            mapQual = "MAPQGE30",
                            norm = "KR"){

  # get all the proper file paths
  # An example path is:
  # K562/100kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_100kb.RAWobserved
  # K562_interchromosomal/100kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_100kb.RAWobserved

  # map resolution in bp to kb/mb substring:
  res2str <- data.frame(
    res=c(1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000),
    resStr=as.character(c("1kb", "5kb", "10kb", "25kb", "50kb", "100kb", "250kb", "500kb", "1mb"))
  )

  # to get proper string representation of the input resolution parameter
  resStr <- as.character(
    res2str$resStr[match(resolution, res2str$res)]
    )
  if(is.na(resStr)){
    stop("Input resolution is not supported. Input resolution: ", resolution,
         " Only the following resolutions are supported: ",
         paste(res2str$res, collapse=", "))
  }
  # example path for GM12878
  # OLD: /project/jgu-cbdm/andradeLab/download/flat/databases/uncompressed/muro/ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GM12878_combined_interchromosomal/50kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_50kb.RAWobserved
  # /project/jgu-cbdm/ibnsalem/translocations/data/Rao2014/GM12878_combined_interchromosomal/10kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_10kb.RAWobserved

  # example for K562:
  # /project/jgu-cbdm/ibnsalem/translocations/data/Rao2014/K562/50kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_50kb.RAWobserved
  # /project/jgu-cbdm/ibnsalem/translocations/data/Rao2014/K562_interchromosomal/50kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_50kb.RAWobserved

  # build path to directory with chromosome subdirectories
  intraDir = file.path(baseDir, ifelse(cell != "GM12878", cell, paste0(cell, "_combined")),
                       paste0(resStr, "_resolution_intrachromosomal"))
  message("INFO: Start parsing from directory: ", intraDir)

  # get all available chromosome names:
  chromosomes = list.dirs(path=intraDir , full.names = FALSE,
                          recursive = FALSE)

  # throw a meaningfull error here if path is wrong
  if( length(chromosomes) == 0 ){
    stop("Could not find chromosomes in ", intraDir)
  }

  # sort
	chromOrder <- match(GenomeInfoDb::seqnames(seqInfo), chromosomes)
	chromosomes <- chromosomes[chromOrder[!is.na(chromOrder)]]

	# create bins for given resolution and chromosome size as GRanges object:
  chromBinGR <- lapply(chromosomes, getBinGR, resolution, seqInfo)
  names(chromBinGR) <- chromosomes

	# combine to binGR for entire genome
	binGR <-unlist(GenomicRanges::GRangesList(chromBinGR))

	# get bin offsets
	binOffset <- getChrToBinOffset(chromosomes, resolution, seqInfo)

	#--------------------------------------------------------------------
	# parse intra-chromosomal interactions
	#--------------------------------------------------------------------

	cisDFlist <- BiocParallel::bplapply(chromosomes, function(chr){

    message(paste("INFO: Begin to parse data for chromosome", chr, "..."))

    # get file path
	  rawInteractionFile = file.path(intraDir, chr, mapQual, paste0(chr, "_", resStr, ".RAWobserved"))

		# parse interaction file
		intData <- data.frame(data.table::fread(rawInteractionFile))

		# get anchor bins as GRanges object
		anchorIdx1 <- binOffset[chr] + (intData[,1] / resolution) +1
		anchorIdx2 <- binOffset[chr] + (intData[,2] / resolution) +1

		# build data.frame with indexes and score
		chrDF <- data.frame(idx1=anchorIdx1, idx2=anchorIdx2, raw=intData[,3])

		return(chrDF)

  })

	#--------------------------------------------------------------------
	# parse inter-chromosomal interactions
	#--------------------------------------------------------------------
  if( !interChromosomal ) {

    intDF <- as.data.frame(data.table::rbindlist(cisDFlist))

  }else{

  	# build path to directory with chromosome subdirectories
  	interDir = file.path(baseDir,
  	                     paste0(
  	                       ifelse(cell != "GM12878", cell, paste0(cell, "_combined")),
  	                       "_interchromosomal"
  	                       ),
  	                     paste0(resStr, "_resolution_interchromosomal"))

  	message("INFO: Start parsing from directory: ", interDir)


  	# get all chromosome pair combination (ordered)
  	chromPairs <- t(utils::combn(chromosomes, 2))

  	transDFlist <- BiocParallel::bpmapply(function(chr1, chr2){

  	  message(paste("INFO: Begin to parse interactions between chromosome", chr1, "and", chr2 , "..."))

  		# get file path
  		# Example file:	#K562_interchromosomal/100kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_100kb.RAWobserved
  	  # baseDir,
  	  # paste0(
  	  #   ifelse(cell != "GM12878", cell, paste0(cell, "_combined")),
  	  #   "_interchromosomal"
  	  # ),
  	  # paste0(resStr, "_resolution_interchromosomal")

      rawInteractionFile = file.path(
        interDir,
        paste0(chr1, "_", chr2),
        mapQual,
        paste0(
          chr1,
          "_",
          gsub("chr", "", chr2),
          "_",
          resStr,
          ".RAWobserved")
        )

  		# parse interaction file
  		intData <- data.frame(data.table::fread(rawInteractionFile))

  		# get anchor bins as GRanges object
  		anchorIdx1 <- binOffset[chr1] + (intData[,1] / resolution) +1
  		anchorIdx2 <- binOffset[chr2] + (intData[,2] / resolution) +1

  		# build data.frame with indexes and score
  		chrPairDF <- data.frame(idx1=anchorIdx1, idx2=anchorIdx2, raw=intData[,3])

  		return(chrPairDF)

  	}, chromPairs[,1], chromPairs[,2], SIMPLIFY=FALSE)

  	# combine all data.frames into a single one
  	intDF <- as.data.frame(data.table::rbindlist(c(cisDFlist, transDFlist)))

  }

	#--------------------------------------------------------------------
	# build GI object
	#--------------------------------------------------------------------
	GI <- InteractionSet::GInteractions(
	  intDF[, 1],
	  intDF[, 2],
	  binGR,
	  raw=intDF[, 3],
	  mode="strict")

	# make user aware of object size
	message(paste(
	  "INFO: Finished parsing of interactions. GI object has size:",
	  format(utils::object.size(GI), units = "auto")
	  ))

  return(GI)
}
