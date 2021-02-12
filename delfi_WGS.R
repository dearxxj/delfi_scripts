###################################################################
# DELFI - DNA EvaLuation of Fragments for early Interception      #
# Adapted from : https://github.com/cancer-genomics/delfi_scripts #
# Date - Feb 09 2021                                              #
# - Input: Bam files from low-coverage (1-2X) WGS                 #
# - Output: Prediction on whether input WGS comes from tumor      #
###################################################################


###########################################
# 0. Parse Options                        #
###########################################

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)

usage <- function() {
	message('Usage: Rscript --vanilla delfi_WGS.R --bam <bam> --out <out>')
}
if (length(options.names) != 2) {
	usage()
	stop('Required options missing!')
}
for (op in options.names) {
	if (!(op %in% c('bam', 'out'))) {
		usage()
		stop(paste0('Unrecognized option --', op))
	}
}

bamfile <- options.args['bam']
outprefix <- options.args['out']


library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(tidyverse)
library(httr)
library(GenomicAlignments)
library(GenomicRanges)
library(getopt)



############################################
## 1. Get hg19 gaps & blacklisted regions  #
############################################
#
#genome <- "hg19"
#mySession <- browserSession()
#genome(mySession) <- genome
#gaps <- getTable(ucscTableQuery(mySession, track="gap"))
#gaps.hg19 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
#					 gaps$chromEnd), type=gaps$type)
#gaps.hg19 <- keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")),
#                           pruning.mode="coarse")
#hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#seqinfo(gaps.hg19) <- seqinfo(hsapiens)[seqlevels(gaps.hg19),]
#
#blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
#blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
#							col_names=c("seqnames", "start",
#							"end", "name", "score"))
#blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
#filters.hg19 <- makeGRangesFromDataFrame(blacklisted.tib,
#                                           keep.extra.columns=TRUE)
#filters.hg19 <- keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")),
#							   pruning.mode="coarse")
#seqinfo(filters.hg19) <- seqinfo(Hsapiens)[seqlevels(filters.hg19),]


###########################################
# 2. Read GAlignmentPairs                 #
###########################################
indexed.bam <- gsub("$", ".bai", bamfile)
if (!file.exists(indexed.bam)) {
    indexBam(bamfile)
}

param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                         isSecondaryAlignment = FALSE,
                                         isUnmappedQuery = FALSE),
                      mapqFilter = 30)
galpdir <- '../'
galp.file <- file.path(galpdir, paste0(outprefix, ".rds"))
galp <- readGAlignmentPairs(bamfile, param = param)
saveRDS(galp, galp.file)

