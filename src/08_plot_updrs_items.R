library(data.table)
library(tidyverse)
library(ggpubr)

setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/violin/eachQ_updrs")
dat <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/df_PD3omics.tsv", sep = '\t') %>%
  select(PATNO, cluster)

pheno <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_longitudinal.tsv", sep = '\t')
pheno_BL <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_BL.tsv", sep = '\t')

# extract PD patients and healthy controls
dat.sub <- pheno_BL %>%
  filter(paste0("PP-", PATNO) %in% dat$PATNO | COHORT == "Healthy Control") %>%
  select(PATNO, COHORT)
# get each Q of UPDRS
pheno2 <- pheno %>%
  filter(PATNO %in% dat.sub$PATNO) %>%
  filter(EVENT_ID == "BL")

merge_df <- merge(dat.sub, pheno2, by = "PATNO") %>%
  mutate(PATNO = paste0("PP-", PATNO))
final_df <- merge(merge_df, dat, by = "PATNO", all.x = T) %>%
  mutate(cluster = ifelse(is.na(cluster), 0, cluster))

n <- table(final_df$cluster)

# library(pals)
# brewer.dark2(8)
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"

my_comparisons <- list(c("0", "1"), c("0", "2"), c("1", "2"))

violin_plot <- function(variable){
  tiff(paste0("violin_", variable, ".tiff"), res = 300, width = 1500, height = 1500)
  p <- ggviolin(final_df, x = "cluster", y = variable,
                color = "cluster",
                add = "boxplot")+
    ylab("UPDRS-I: Daytime sleepiness")+
    scale_color_brewer(palette = "Dark2")+
    scale_x_discrete(labels = c("0" = paste0("HC\n (n = ", n[1], ")"), 
                                "1" = paste0("C1\n (n = ", n[2], ")"), 
                                "2" = paste0("C2\n (n = ", n[3], ")")))+
    stat_compare_means(comparisons = my_comparisons, size = 4)+
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20))
  print(p)
  dev.off()
}

for (i in 8:74){
  outcome <- colnames(final_df)
  violin_plot(variable = outcome[i])
}




