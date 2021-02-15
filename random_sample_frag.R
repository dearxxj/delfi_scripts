# Random sample fragment counts 100 time for tumor and control
# to generate 100 tumor fragment profile and 100 healthy fragment profile
require(tidyverse)
require(GenomicRanges)

AB_cancer <- readRDS('result/K562_bin_100kb.rds')
AB_healthy <- readRDS('result/NA12878_bin_100kb.rds')

outdir <- "../random_100kb"
dir.create(outdir, showWarnings=F)

sample_ref <- read_csv("sample_reference.csv")

for (i in 1:nrow(sample_ref)) {
	if (sample_ref[i, 'Patient Type'] != "Healthy") {
		AB_random <- AB_cancer[, 1:18]
	} else {
		AB_random <- AB_healthy[, 1:18]
	}
	random_frags <- AB_random %>% as_tibble %>% select(short, long, ratio, nfrags, short.corrected, long.corrected, nfrags.corrected, ratio.corrected, frag.gc) %>% sample_n(length(AB_random))
	AB_random$short <- random_frags$short
	AB_random$long <- random_frags$long
	AB_random$ratio <- random_frags$ratio
	AB_random$nfrags <- random_frags$nfrags
	AB_random$short.corrected <- random_frags$short.corrected
	AB_random$long.corrected <- random_frags$long.corrected
	AB_random$ratio.corrected <- random_frags$ratio.corrected
	AB_random$nfrags.corrected <- random_frags$nfrags.corrected
	AB_random$frag.gc <- random_frags$frag.gc
	
	set.seed(i)
	id <- sample_ref[i, 'WGS ID']
	saveRDS(AB_random, file.path(outdir, paste0(id, '_bin_100kb.rds')))
}
