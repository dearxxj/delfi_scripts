# Random sample fragment counts 100 time for tumor and control
# to generate 100 tumor fragment profile and 100 healthy fragment profile
require(dplyr)

AB_cancer <- load('../K562_bin_100kb.rds')
AB_healthy <- load('../NA12878_bin_100kb.rds')
AB_cancer <- AB_cancer %>% select(-matches("X"))
AB_healthy <- AB_healthy %>% select(-matches("X"))

outdir <- "../random_100kb"

sample_ref <- read_csv("sample_reference.csv")
sample_cancer <- sample_ref %>% filter(`Patient Type` != "Healthy") %>% head(100)
sample_healthy <- sample_ref %>% filter(`Patient Type` == "Healthy") %>% head(100)
sample_select <- bind_rows(sample_cancer, sample_healthy)

write_csv(sample_select, 'sample_select200.csv')


