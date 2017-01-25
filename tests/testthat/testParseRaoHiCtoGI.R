context("Parse Hi-C of Rao et al 2014")

# example directory structure in extdata
baseDir <- system.file("extdata", "Rao2014", package="chromint")

test_that("Binning of single chromosome has correct number of bins and no overlap", {

  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  binGR <- getBinGR("chr1", 50000, seqinfo(txdb))

  # get seqlength from seqinfo object
  chrLen <- seqlengths(seqinfo(txdb))[["chr1"]]

  expect_equal(ceiling(chrLen/50000), length(binGR))
  expect_equal(length(findOverlaps(binGR, binGR)), length(binGR))
})


test_that("bin offset of cromosomes is increasing", {

  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  offset <- getChrToBinOffset(c("chr1", "chr2"), 50000, seqinfo(txdb))

  expect_equal(offset[["chr1"]], 0)

  curOff <- 0
  for(i in 1:length(offset)){
    expect_equal(offset[[i]] >= curOff, TRUE)
    curOff <- offset[[i]]
  }
})


test_that("Parsing for K562 cell gives interactions", {

  seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

  cell <- "K562"
  resolution <- 50000

  gi <- parseRaoHiCtoGI(cell, resolution, baseDir, seqInfo)

  expect_equal(length(gi) > 0, TRUE)
  expect_equal(max(width(InteractionSet::regions(gi))), resolution)
})

test_that("Parsing for special case of GM12878 cell gives interactions", {

  seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)

  cell <- "GM12878"
  resolution <- 50000

  gi <- parseRaoHiCtoGI(cell, resolution, baseDir, seqInfo)

  expect_equal(length(gi) > 0, TRUE)
  expect_equal(max(width(InteractionSet::regions(gi))), resolution)
})
