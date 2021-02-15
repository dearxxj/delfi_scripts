###################################################################
# DELFI - DNA EvaLuation of Fragments for early Interception      #
# Adapted from : https://github.com/cancer-genomics/delfi_scripts #
# Date - Feb 09 2021                                              #
# Params:                                                         #
#   --bam: A Bam file from low-coverage (1-2X) WGS                #
#   --out: output prefix                                          #
# Output: A fragment profile at 100kb bins for the bam            #
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
	message('Usage: Rscript --vanilla delfi_WGS.R --bam <bam> --out_prefix <output_prefix> --out_dir <output_directory>')
}
if (length(options.names) != 2) {
	usage()
	stop('Required options missing!')
}
for (op in options.names) {
	if (!(op %in% c('bam', 'out_prefix', 'out_dir'))) {
		usage()
		stop(paste0('Unrecognized option --', op))
	}
}

bamfile <- options.args['bam']
outprefix <- options.args['out_prefix']
outdir <- options.args['out_dir']


library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(tidyverse)
library(httr)
library(GenomicAlignments)
library(GenomicRanges)
library(getopt)
library(Homo.sapiens)
library(biovizBase)
library(RCurl)
library(readxl)
library(multidplyr)



###########################################
# 1. Get hg19 gaps & blacklisted regions  #
###########################################

genome <- "hg19"
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
gaps.hg19 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
					 gaps$chromEnd), type=gaps$type)
gaps.hg19 <- keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")),
                           pruning.mode="coarse")
hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
seqinfo(gaps.hg19) <- seqinfo(hsapiens)[seqlevels(gaps.hg19),]
save(gaps.hg19, file=file.path(outdir, "gaps.hg19.rda"))

blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
							col_names=c("seqnames", "start",
							"end", "name", "score"))
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
filters.hg19 <- makeGRangesFromDataFrame(blacklisted.tib,
                                           keep.extra.columns=TRUE)
filters.hg19 <- keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")),
							   pruning.mode="coarse")
seqinfo(filters.hg19) <- seqinfo(Hsapiens)[seqlevels(filters.hg19),]
save(filters.hg19, file=file.path(outdir, "filters.hg19.rda"))


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
galp.file <- file.path(outdir, paste0(outprefix, ".rds"))
galp <- readGAlignmentPairs(bamfile, param = param)
saveRDS(galp, galp.file)


###########################################
# 3. Filter fragments / Get GC content    #
###########################################
galp <- readRDS(galp.file)
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
             on.discordant.seqnames="drop")

## filter outliers
w.all <- width(frags)
q.all <- quantile(w.all, c(0.001, 0.999))
frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]

frags$gc <- GCcontent(Hsapiens, unstrand(frags))

frag.file <- file.path(outdir, paste0(outprefix, "_frags.rds"))
saveRDS(frags, frag.file)


###########################################
# 4. Bin compartments, GC correction      #
###########################################

gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage, na.rm=T)
}


filename <- file.path(outdir, paste0(outprefix, "_bin_100kb.rds"))

load("./filters.hg19.rda")
load("./gaps.hg19.rda")

ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

AB <- read.table(textConnection(ABurl), header=TRUE)
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))

tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]

arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)

## These bins had no coverage
#AB <- AB[-c(8780, 13665)]
fragments <- readRDS(frag.file)

### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg19))]
w.all <- width(fragments)

fragments <- fragments[which(w.all >= 100 & w.all <= 220)]
w <- width(fragments)

frag.list <- split(fragments, w)

counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
if(min(w) > 100) {
    m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
                 dimnames=list(rownames(counts), 100:(min(w)-1)))
    counts <- cbind(m0, counts)
}

olaps <- findOverlaps(fragments, AB)
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
bingc <- rep(NA, length(AB))
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))

### Get modes
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

short <- rowSums(counts[,1:51])
long <- rowSums(counts[,52:121])
ratio <- short/long
ratio[is.na(ratio)] <- NA
ratio[is.infinite(ratio)] <- NA
short.corrected=gc.correct(short, bingc)
long.corrected=gc.correct(long, bingc)
nfrags.corrected=gc.correct(short+long, bingc)
ratio.corrected=gc.correct(ratio, bingc)

AB$short <- short
AB$long <- long
AB$ratio <- ratio
AB$nfrags <- short+long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc

for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]

saveRDS(AB, filename)
