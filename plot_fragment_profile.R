## Fragment profile plot
df <- readRDS("result/bins_5mbcompartments.rds")
df <- df %>% dplyr::filter(type %in% c('Healthy', 'Lung'))
health_sample <- df %>% filter(type=="Healthy") %>% dplyr::select(sample) %>% unique %>% head(12)
lung_sample <- df %>% filter(type=="Lung") %>% dplyr::select(sample) %>% unique %>% head(12)
df <- df %>% dplyr::filter(sample %in% c(health_sample$sample, lung_sample$sample))

df.fr <- df %>% group_by(sample) %>% mutate(ratio.centered = scale(ratio.corrected2, scale=F))

df.fr %>% group_by(type) %>% filter(bin==1) %>% summarize(n=n())
df.fr$type <- relevel(factor(df.fr$type), "Healthy")
tissue <- c(Healthy = "NA12878 WGS", Lung = "K562 WGS")

arm <- df.fr %>% group_by(arm) %>% summarize(n=n()) %>% mutate(arm = as.character(arm))
small.arms <- setNames(c("", "12q", "", "16q", "", "17q", "", "18q", "", "19", "", "20", "21", "22"), c("12p", "12q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

mytheme <- theme_classic(base_size=12) + theme( axis.text.x = element_blank(),
                              axis.ticks.x=element_blank(),
                           strip.background = element_blank(),
                           strip.text.y = element_text(face="bold", size=50),
                           strip.text.x = element_text(face="bold", size=30),
                           axis.title.x = element_text(face="bold", size=30),
                           axis.title.y = element_text(face="bold", size=50),
                           axis.text.y = element_text(size = 20),
                           plot.title = element_text(size=15),
						   legend.position = "none",
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()	)
gg <- ggplot(df.fr, aes(x=bin, y=ratio.centered, group=sample)) + geom_line(color="gray50", size=1)
gg <- gg + labs(x="", y="Fragmentation profile", color="", size=1)
gg <- gg + facet_grid(type~arm, switch="both",space="free_x", scales="free_x", labeller=labeller(type=tissue, arm=arm.labels))
gg <- gg + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
gg <- gg  + mytheme
gg
save_plot("result/fragment_profile.png", gg, ncol=1, nrow=1, base_height=12, base_width=35)


## ROC curvy for prediction
prediction <- read_csv("result/predictions_gbm.csv")
y <- prediction$type
y[y != "Healthy"] <- "Cancer"
par(pty='s')
png("result/prediction_ROC.png", height = 480, width = 480)
plot(roc(y, prediction$Cancer), print.auc=T, print.auc.x=0.8,print.auc.y=0.8,legacy.axes=T)
dev.off()
