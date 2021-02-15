###################################################################
# DELFI - DNA EvaLuation of Fragments for early Interception      #
# Adapted from : https://github.com/cancer-genomics/delfi_scripts #
# Date - Feb 09 2021                                              #
# Params:                                                         #
#    --dir: directory of 100kb binned fragment profile            #
#    --sample: sample reference csv file                          #
#    --out: output prefix                                         #
# Output:             #
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
	message('Usage: Rscript --vanilla summarize_data_gbm.R --bindir <100kb bin directory> --sample <sample_reference.csv> --outdir <output directory>')
}
if (length(options.names) != 3) {
	usage()
	stop('Required options missing!')
}
for (op in options.names) {
	if (!(op %in% c('bindir', 'sample', 'outdir'))) {
		usage()
		stop(paste0('Unrecognized option --', op))
	}
}

bindir <- options.args['bindir']
sample_ref <- options.args['sample']
outdir <- options.args['outdir']

library(tidyverse)
library(GenomicRanges)
library(multidplyr)
library(readxl)
library(caret)
library(gbm)
library(pROC)

###############################################
# 5. Aggregate frag profile for all samples   #
###############################################

metadata <- read_csv(sample_ref)
ids <- metadata %>% select(`WGS ID`) %>% unlist()
files <- file.path(bindir, paste0(ids, "_bin_100kb.rds"))

bins.list <- lapply(files, readRDS)
tib.list <- lapply(bins.list, as_tibble)
names(tib.list) <- ids
tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%
    bind_rows() %>% select(id, everything())

tib.list <- tib.list %>% select(-matches("X"))
saveRDS(tib.list, file.path(outdir, "bins_100kbcompartments.rds"))

df.fr <- readRDS(file.path(outdir, "bins_100kbcompartments.rds"))
master <- read_csv(sample_ref)

df.fr2 <- inner_join(df.fr, master, by=c("id"="WGS ID"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)

## combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric end and remove the bin closest to the centromere if it is
## smaller than 5mb.
df.fr2 <- df.fr2 %>% group_by(id, arm) %>%
    mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                           ceiling(rev((1:length(arm))/50) )))

df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
    summarize(short2=sum(short, na.rm=T),
              long2=sum(long, na.rm=T),
              short.corrected2=sum(short.corrected, na.rm=T),
              long.corrected2=sum(long.corrected, na.rm=T),
              hic.eigen=mean(eigen, na.rm=T),
              gc=mean(C.G, na.rm=T),
              ratio2=mean(ratio, na.rm=T),
              ratio.corrected2=mean(ratio.corrected, na.rm=T),
              nfrags2=sum(nfrags, na.rm=T),
              nfrags.corrected2=sum(nfrags.corrected, na.rm=T),
              domain = median(as.integer(domain), na.rm=T),
              short.var=var(short.corrected, na.rm=T),
              long.var=var(long.corrected, na.rm=T),
              nfrags.var=var(nfrags.corrected, na.rm=T),
              mode_size=unique(mode),
              mean_size=unique(mean),
              median_size=unique(median),
              q25_size=unique(quantile.25),
              q75_size=unique(quantile.75),
              start=start[1],
              end=rev(end)[1],
              binsize = n()) %>% dplyr::rename("sample" = "id")
### assign bins
df.fr3 <- inner_join(df.fr3, master, by=c("sample"="WGS ID"))
df.fr3 <- df.fr3 %>% mutate(type = gsub(" Cancer|carcinoma", "", `Patient Type`, ignore.case=TRUE))

df.fr3 <- df.fr3 %>% filter(binsize==50)
df.fr3 <- df.fr3 %>% group_by(sample) %>% mutate(bin = 1:length(sample))

saveRDS(df.fr3, file.path(outdir, "bins_5mbcompartments.rds"))

###############################################
# 6. Summarize data and prepare input for gbm #
###############################################
master <- read_csv(sample_ref)
df.fr3 <- readRDS(file.path(outdir, "bins_5mbcompartments.rds"))

healthy.median <- df.fr3 %>%
    group_by(bin) %>% 
    summarize(median.cov=median(nfrags2, na.rm=TRUE),
              median.short=median(short2, na.rm=TRUE),
              median.long=median(long2, na.rm=TRUE),
              median.ratio=median(ratio2, na.rm=TRUE),
              median.corrected.cov=median(nfrags.corrected2, na.rm=TRUE),
              median.corrected.short=median(short.corrected2, na.rm=TRUE),
              median.corrected.long=median(long.corrected2, na.rm=TRUE),
              median.corrected.ratio=median(ratio.corrected2, na.rm=TRUE),
              median.corrected.ratio2=median(short.corrected2/long.corrected2, na.rm=TRUE))
summary.df <- df.fr3 %>% ungroup() %>% group_by(sample, type) %>%
    summarize(cov.cor=cor(nfrags2, healthy.median$median.cov, method="pearson", use="complete.obs"),
              short.cor=cor(short2, healthy.median$median.short, method="pearson", use="complete.obs"),
              long.cor=cor(long2, healthy.median$median.long, method="pearson", use="complete.obs"),
              ratio.cor=cor(ratio2, healthy.median$median.ratio, method="pearson", use="complete.obs"),
              cov.corrected.cor=cor(nfrags.corrected2, healthy.median$median.corrected.cov, method="pearson", use="complete.obs"),
              short.corrected.cor=cor(short.corrected2, healthy.median$median.corrected.short, method="pearson", use="complete.obs"),
              long.corrected.cor=cor(long.corrected2, healthy.median$median.corrected.long, method="pearson", use="complete.obs"),
              ratio.corrected.cor=cor(ratio.corrected2, healthy.median$median.corrected.ratio, method="pearson", use="complete.obs"),
              ratio2.corrected.cor=cor(short.corrected2/long.corrected2, healthy.median$median.corrected.ratio2, method="pearson", use="complete.obs"),
              nfrags = sum(nfrags2),
              mode_size=unique(mode_size),
              mean_size=unique(mean_size),
              median_size=unique(median_size),
              q25_size=unique(q25_size),
              q75_size=unique(q75_size),
              hqbases_analyzed = 100*sum(nfrags)*2,
              coverage = hqbases_analyzed/(504*5e6)
              )

summary.df <- inner_join(summary.df, master, by=c("sample"="WGS ID"))
summary.df$`type` = relevel(as.factor(summary.df$`type`), "Healthy")

saveRDS(summary.df, file.path(outdir, "summary_tibble.rds"))


##################################################################
# 7. Prediction of tumore/healthy with Stochastic Gradient Boost #
##################################################################


df.fr3 <- readRDS(file.path(outdir, "bins_5mbcompartments.rds"))
summary.df <- readRDS(file.path(outdir, "summary_tibble.rds"))

features.cov <- df.fr3  %>% ungroup() %>%
    select(nfrags.corrected2, sample, bin) %>%
    spread(sample, nfrags.corrected2) %>%
    select(-bin) %>% 
    na.omit() %>%
    scale() %>%
    t() %>%
    as.data.frame()

features.short <- df.fr3  %>% ungroup() %>%
    select(short.corrected2, sample, bin) %>%
    spread(sample, short.corrected2) %>%
    select(-bin) %>% 
    na.omit() %>%
    scale() %>%
    t() %>%
    as.data.frame()

features.sl <- cbind(features.cov, features.short)
colnames(features.sl) <- c(paste0("total", 1:ncol(features.cov)), paste0("short", 1:ncol(features.short)))
features.sl$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

features <- cbind(features.sl,
             as.matrix(summary.df %>% ungroup() %>%
                       select(contains("Z Score"))))
features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
#                      preProcOptions=list(thres = 0.90),
                     summaryFunction = twoClassSummary)
set.seed(1234)
model_gbm <- caret::train(type ~ .,
                               data = features,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150,
                                                   interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

####### Only short/total coverage
set.seed(1234)
model_sl <- caret::train(type ~ .,
                               data = features.sl,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

###### Only z-scores
features.z <- summary.df %>% ungroup() %>% select(contains("Z Score"))
features.z$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

set.seed(1234)
model_z <- caret::train(type ~ .,
                               data = features.z,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150,
                                                   interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

#### Save
models.list <- list("all"=model_gbm, "SL"=model_sl, "z"=model_z)
saveRDS(models.list, file.path(outdir, "models_list.rds"))

pred.tbl <- model_gbm$pred %>% filter(n.trees==150, interaction.depth==3) %>%
group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample <- rownames(features)
pred.tbl <- inner_join(pred.tbl, summary.df)

## 95% specificity
cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
           arrange(desc(Cancer)))$Cancer[11]
cutoff98 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[5]
## 90% specificity cutoff to be used in tissue prediction.
cutoff90 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[21]

pred.tbl <- pred.tbl %>%
    mutate(detected95 = ifelse(Cancer > cutoff, "Detected", "Not detected"),
       detected98 = ifelse(Cancer > cutoff98, "Detected", "Not detected"),
       detected90 = ifelse(Cancer > cutoff90, "Detected", "Not detected"),
       stage = gsub("A|B|C", "", `Stage at Diagnosis`))

write.csv(inner_join(summary.df %>% select(-contains("Z Score")), pred.tbl %>%
                     select(rowIndex, sample, stage, Cancer, detected95, detected98),
                     by=c("sample"="sample")), file.path(outdir, "predictions_gbm.csv"),
          row.names=FALSE)
