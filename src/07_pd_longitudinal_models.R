library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(ggeffects)

df_pheno <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/df_PD3omics.tsv", sep = '\t')
longit_pheno <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_longitudinal.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>%
  arrange(PATNO)

df <- merge(longit_pheno, df_pheno, by = "PATNO") %>%
  filter(time %in% c(0, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48, 54, 60))

##### linear mixed model #####
df$type_mutation <- factor(df$type_mutation, levels = c("nomut", "GBA", "LRRK2"))
lmm <- lmer(MCATOT.x ~ SEX + AGE_AT_VISIT.y + duration_yrs + type_mutation + time*cluster + (1|PATNO), data = df)
summary(lmm)
pred.lmm <- ggpredict(lmm, terms = c("time", "cluster"))

#' Adaptive palette (discrete).
#' Create a discrete palette which will use the first n colors from
#' the supplied color values, and interpolate after n.
adaptive_pal <- function(values) {
  force(values)
  function(n = 10) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

pal_npg_adaptive <- function(palette = c("nrc"), alpha = 1) {
  palette <- match.arg(palette)
  
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db$"npg"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  adaptive_pal(unname(alpha_cols))
}

scale_color_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("colour", "npg", pal_npg_adaptive(palette, alpha), ...)
}

scale_fill_npg_adaptive <- function(palette = c("nrc"), alpha = 1, ...) {
  palette <- match.arg(palette)
  discrete_scale("fill", "npg", pal_npg_adaptive(palette, alpha), ...)
}

df$group <- as.factor(df$cluster)
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/longitudinal/GI.tiff", res = 300, width = 2500, height = 1500)
ggplot(pred.lmm, aes(x, predicted, colour = group))+
  geom_line()+
  # facet_wrap(~ group)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.3, colour = NA)+
  scale_color_manual(values = c("1"="coral2",
                                "2"="darkcyan"))+
  scale_fill_manual(values = c("1"="coral2",
                               "2"="darkcyan"))+
  geom_point(data = df, aes(x = time, y = GI.x, colour = group, alpha = 0.3))+
  # scale_color_npg_adaptive()+
  # scale_fill_npg_adaptive()+
  xlab("Months")+
  ylab("GI")+
  theme(
    # legend.position = "none", 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15))
dev.off()

