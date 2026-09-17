library(iClusterPlus)
library(data.table)
library(tidyverse)
library(ggpubr)
library(lattice)

meta <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/metadata/PPMI_Curated_Data_Cut_Public_20240129.tsv") %>%
  filter(EVENT_ID == "BL" & type_mutation == "GBA1") %>%
  mutate(log2_ttau = log2(tau + 1),
         log2_ptau = log2(ptau + 1),
         log2_asyn = log2(asyn + 1),
         log2_abeta = log2(abeta + 1)) %>% 
  arrange(PATNO)


# CSF Proteomics
setwd('/Volumes/Extreme SSD/work_bidmc/PPMI/proteomics/p151/Proteomic_Analysis')
anno <- read.csv('PPMI_Project_151_pqtl_Analysis_-Table 1.csv')
proteomics_df <- fread('proteomics_220621.tsv', sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>% arrange(PATNO) %>%
  column_to_rownames(var = "PATNO") %>%
  as.data.frame()
probe_human <- anno %>%
  filter(ORGANISM == "Human") %>%
  select(SOMA_SEQ_ID)
proteomics_human <- proteomics_df[, colnames(proteomics_df) %in% probe_human$SOMA_SEQ_ID] #1158x4779
name_symbol <- anno[match(colnames(proteomics_human), anno[, 'SOMA_SEQ_ID']), 'TARGET_GENE_SYMBOL'] #4779 SOMAmers targeting 4102 unique proteins
colnames(proteomics_human) <- name_symbol

# Urine Proteomics
urine <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/proteomics/Urine_Proteomics/DF2_log2_pc4.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>% arrange(PATNO) %>%
  column_to_rownames(var = "PATNO") %>%
  as.data.frame() %>%
  select(-c(1:9))

## QC: exclude proteins with > 20% missing values
NAcount_col <- function(df){
  apply(df, 2, function(x){ sum(is.na(x)) })
}

urine_qc <- urine %>% 
  select(-which(NAcount_col(.) > nrow(.)*0.2)) # 979x2170
## impute missing values with the mean value
for(j in 1:ncol(urine_qc)){
  urine_qc[is.na(urine_qc[, j]), j] <- mean(urine_qc[, j], na.rm = T)
}

# Metabolics
dat5 <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/metabolomics_PPMI180/Metabolomic_Analysis/Metabolomic_Analysis_of_LRRK2_PD_5_of_5.csv") %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>%
  arrange(PATNO, CLINICAL_EVENT) %>%
  filter(UNITS == 'adjusted_area_ratio') %>%
  select(PATNO, SEX, COHORT, CLINICAL_EVENT, TESTNAME, TESTVALUE)
dat5_wide <- dat5 %>%
  pivot_wider(names_from = TESTNAME, values_from = TESTVALUE) %>%
  distinct(PATNO, .keep_all = T) # 607x302
metabolics <- dat5_wide %>%
  column_to_rownames(var = "PATNO") %>% select(-c(1:4)) # 607x297

# QC: exclude metabolites with > 20% missing values
metabolics_qc <- metabolics %>% select(-which(NAcount_col(.) > nrow(.)*0.2)) # 607x280
## impute missing value with the mean value
for(j in 1:ncol(metabolics_qc)){
  metabolics_qc[is.na(metabolics_qc[, j]), j] <- mean(metabolics_qc[, j], na.rm = T)
}
metabolics_qc <- log2(metabolics_qc + 1)

# Transcriptomics
rna <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/rnaseq/log2CPM_top2000_chr22_coding.tsv") %>%
  pivot_longer(-Gene) %>%
  pivot_wider(names_from = Gene, values_from = value) %>%
  arrange(name) %>%
  column_to_rownames(var = "name")

comm_id <- intersect(meta$PATNO, rownames(proteomics_human)) #225
comm_id <- intersect(comm_id, rownames(rna)) #116
comm_id <- intersect(comm_id, rownames(urine_qc)) #109
# comm_id <- intersect(comm_id, rownames(metabolics_qc)) #47

df.csf_proteomics <- as.matrix(proteomics_human[rownames(proteomics_human) %in% comm_id,])
df.transcriptomics <- as.matrix(rna[rownames(rna) %in% comm_id,])
df.urine_proteomics <- as.matrix(urine_qc[rownames(urine_qc) %in% comm_id,])
# df.metabolics <- as.matrix(metabolics_qc[rownames(metabolics_qc) %in% comm_id,])

# most variable features
compute_cv <- function(x){ sd(x)/mean(x) }
cv.csf_proteomics <- apply(df.csf_proteomics, 2, compute_cv)
cv.csf_proteomics <- sort(cv.csf_proteomics, decreasing = T)
cv.csf_proteomics_top2000 <- names(head(cv.csf_proteomics, 2000))
cv.csf_proteomics_top2000 <- df.csf_proteomics[, cv.csf_proteomics_top2000]
colnames(cv.csf_proteomics_top2000) <- make.unique(colnames(cv.csf_proteomics_top2000), sep = '.')

cv.urine_proteomics <- apply(df.urine_proteomics, 2, compute_cv)
cv.urine_proteomics <- sort(cv.urine_proteomics, decreasing = T)
cv.urine_proteomics_top2000 <- names(head(cv.urine_proteomics, 2000))
cv.urine_proteomics_top2000 <- df.urine_proteomics[, cv.urine_proteomics_top2000]

################################################################################
#                                iClusterBayes                                 #
################################################################################
set.seed(123)
date()
bayfit <- tune.iClusterBayes(cpus = 6, 
                             # dt1 = df.metabolics,
                             dt1 = cv.csf_proteomics_top2000,
                             dt2 = cv.urine_proteomics_top2000,
                             dt3 = df.transcriptomics,
                             type = c('gaussian', 'gaussian', 'gaussian'),
                             K = 1:6, n.burnin = 12000, # number of MCMC burnin
                             n.draw = 22000, # number of MCMC draw
                             prior.gamma = rep(0.1, 3),
                             sdev = 0.5, pp.cutoff = 0.5)
date()
save.image(file = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/gba1_3omics.RData")
# load("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/gba1_3omics.RData")

################################################################################
#             Model selection using Bayesian information criteria (BIC)        #
################################################################################
allBIC = NULL
devratio = NULL
nK = length(bayfit$fit)
for(i in 1:nK){
  allBIC = c(allBIC, bayfit$fit[[i]]$BIC)
  devratio = c(devratio,bayfit$fit[[i]]$dev.ratio)
}

# plot the number of clusters vs. percent of explained variation
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/BIC_deviance_gba1_3omics.tiff")
par(mar=c(4.0, 4.0, 0.5, 0.5), mfrow=c(1, 2))
plot(1:nK, allBIC, type="b", xlab="k", ylab="BIC", pch=c(19, 1, 1, 1, 1, 1))
plot(1:nK, devratio, type="b", xlab="k", ylab="Deviance ratio", pch=c(19, 1, 1, 1, 1, 1))
dev.off()

k <- which.max(devratio)
best.cluster.Bayes = bayfit$fit[[k]]$clusters
table(best.cluster.Bayes)

df.pheno <- meta[meta$PATNO %in% comm_id,] %>%
  mutate(cluster = factor(best.cluster.Bayes))

tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/posterior_probability_3omics.tiff", res = 300, height = 2500, width = 1500)
par(mfrow = c(3, 1))
# plot(bayfit$fit[[k]]$beta.pp[[1]], xlab = "Feature", ylab = "Posterior probability", main = "Plasma Metabolomics")
plot(bayfit$fit[[k]]$beta.pp[[1]], xlab = "Feature", ylab = "Posterior probability", main = "CSF Proteomics")
plot(bayfit$fit[[k]]$beta.pp[[2]], xlab = "Feature", ylab = "Posterior probability", main = "Urine Proteomics")
plot(bayfit$fit[[k]]$beta.pp[[3]], xlab = "Feature", ylab = "Posterior probability", main = "WB Transcriptomics")
dev.off()

features <- alist()
# features[[1]] <- colnames(df.metabolics)
features[[1]] <- colnames(cv.csf_proteomics_top2000)
features[[2]] <- colnames(cv.urine_proteomics_top2000)
features[[3]] <- colnames(df.transcriptomics)

signatures <- alist()
for(i in 1:3){
  signatures[[i]] <- (features[[i]])[which(bayfit$fit[[k]]$beta.pp[[i]] > 0.5)]
}
# signatures[[1]] <- (features[[1]])[which(bayfit$fit[[k]]$beta.pp[[1]] > 0.5)]
# signatures[[2]] <- (features[[2]])[which(bayfit$fit[[k]]$beta.pp[[2]] > 0.99)]
# signatures[[3]] <- (features[[3]])[which(bayfit$fit[[k]]$beta.pp[[3]] > 0.5)]
# signatures[[4]] <- (features[[4]])[which(bayfit$fit[[k]]$beta.pp[[4]] > 0.5)]
length(signatures[[1]])
length(signatures[[2]])
length(signatures[[3]])
# length(signatures[[4]])

# ### venn diagram ###
# library(VennDiagram)
# csf <- signatures[[1]]
# # urine <- signatures[[3]]
# # urine <- gsub("\\_.*|\\;.*", "", urine)
# rna <- signatures[[2]]
# 
# v <- venn.diagram(x = list(csf, urine, rna),
#                   category.names = c('CSF_proteomics', 'Urine_proteomics', 'WB_rna'),
#                   filename = NULL,
#                   cex=1.5,
#                   cat.pos = c(135, 180, 180),
#                   cat.dist = c(0.01, 0.01, 0.025),
#                   margin = 0.05,
#                   col = c("#440154ff", "#31B7BC99", "#F39200FF"),
#                   fill = c("#440154ff", "#31B7BC99", "#F39200FF"))
# # col = c("#440154ff", "#21908dff", "#F39200FF"),
# # fill = c("#440154ff", "#21908dff", "#F39200FF"))
# tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD_4omics/venn_PD_4omics.tiff", res = 300, width = 1500, height = 1500)
# grid.newpage()
# grid.draw(v)
# dev.off()


################################################################################
#                               Generate Heatmap                               #
################################################################################
# Ggenomic features with posterior probability > 0.5 were selected
# truncate the values for better image plot
trunc_func <- function(x, threshold){
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  return(x)
}

# scale.metabolomics <- scale(df.metabolics)
scale.csf_proteomics <- scale(cv.csf_proteomics_top2000)
scale.urine_proteomics <- scale(cv.urine_proteomics_top2000)
scale.transcriptomics <- scale(df.transcriptomics)

# trunc.metabolomics <- trunc_func(scale.metabolomics, 4)
trunc.csf_proteomics <- trunc_func(scale.csf_proteomics, 2)
trunc.urine_proteomics <- trunc_func(scale.urine_proteomics, 5)
trunc.transcriptomics <- trunc_func(scale.transcriptomics, 4)

# col.scheme = alist()
# col.scheme[[1]] = bluered(256)
# col.scheme[[2]] = bluered(256)
# col.scheme[[3]] = bluered(256)
# col.scheme[[4]] = bluered(256)
# tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/lrrk2/heatmap_lrrk2_4omics.tiff", res = 300, height = 1500, width = 2500)
# plotHMBayes(fit = bayfit$fit[[k]],
#             datasets = list(trunc.metabolomics, trunc.csf_proteomics, trunc.urine_proteomics, trunc.transcriptomics),
#             type = c("gaussian", "gaussian", "gaussian", "gaussian"),
#             threshold = c(0.2, 0.5, 0.8, 0.5),
#             col.scheme = col.scheme,
#             row.order = rep(T, 4),
#             scale = rep(F, 4),
#             sparse = rep(T, 4),
#             cap = rep(T, 4))
# dev.off()


# Heatmaps
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gplots)

df.pheno$scale.ttau <- scale(df.pheno$log2_ttau)
df.pheno$scale.ptau <- scale(df.pheno$log2_ptau)
df.pheno$scale.asyn <- scale(df.pheno$log2_asyn)
df.pheno$scale.abeta <- scale(df.pheno$log2_abeta)

df.pheno$SEX <- ifelse(df.pheno$SEX == 1, "Male", "Female")

df.pheno$saa[df.pheno$CSFSAA_new == 0] <- "Negative"
df.pheno$saa[df.pheno$CSFSAA_new == 1] <- "Positive"

df.pheno$PD <- ifelse(df.pheno$cohort == 1, "PD", 
                      ifelse(df.pheno$cohort == 2, "HC", "Prodromal"))

col_markers <- list(cluster = c("1" = "#C7522B", "2" = "#089392"),
                    PD_status = c("HC" = "#24693D", "PD" = "#AE123A", "Prodromal" = "#66ACCF"),
                    # SEX = c("Female" = "#E49183", "Male" = "#66ACCF"),
                    # mutation = c("GBA" = "#F3920099", "LRRK2" = "#D5131799", "SNCA" = "#F4D166", "sPD" = "#1B9E77"),
                    # SAA = c("Positive" = "#C57348", "Negative" = "#63B1C1"),
                    # GI = colorRamp2(range(df.pheno$scopa_gi, na.rm = T), colors = c("#ECD999", "#EF4868")),
                    # ttau = colorRamp2(range(df.pheno$scale.ttau, na.rm = T), colors = c("#F4D166", "#24693D")),
                    # ptau = colorRamp2(range(df.pheno$scale.ptau, na.rm = T), colors = c("#F4D166", "#24693D")),
                    # asyn = colorRamp2(range(df.pheno$scale.asyn, na.rm = T), colors = c("#F4D166", "#24693D")),
                    abeta = colorRamp2(range(df.pheno$scale.abeta, na.rm = T), colors = c("#F4D166", "#24693D"))
                    # updrs1 = colorRamp2(range(df.pheno$updrs1_score, na.rm = T), colors = c("#FFBEB2", "#AE123A")),
                    # updrs2 = colorRamp2(range(df.pheno$updrs2_score, na.rm = T), colors = c("#FFBEB2", "#AE123A")),
                    # updrs3 = colorRamp2(range(df.pheno$updrs3_score, na.rm = T), colors = c("#FFBEB2", "#AE123A")),
                    # MoCA = colorRamp2(range(df.pheno$moca, na.rm = T), colors = c("#FFBEB2", "#AE123A")),
                    # caudate = colorRamp2(range(df.pheno$mean_caudate, na.rm = T), colors = c("#FFC685", "#9E3D22")),
                    # putamen = colorRamp2(range(df.pheno$mean_putamen, na.rm = T), colors = c("#FFC685", "#9E3D22")),
                    # striatum = colorRamp2(range(df.pheno$mean_striatum, na.rm = T), colors = c("#FFC685", "#9E3D22"))
                    )
col_ha <- HeatmapAnnotation(cluster = df.pheno$cluster,
                            PD_status = df.pheno$PD,
                            # SEX = df.pheno$SEX,
                            # SAA = df.pheno$saa,
                            # GI = df.pheno$scopa_gi,
                            # ttau = df.pheno$scale.ttau, 
                            # ptau = df.pheno$scale.ptau,
                            # asyn = df.pheno$scale.asyn, 
                            abeta = df.pheno$scale.abeta,
                            # caudate = df.pheno$mean_caudate, 
                            # putamen = df.pheno$mean_putamen, 
                            # striatum = df.pheno$mean_striatum,
                            # updrs1 = df.pheno$updrs1_score,
                            # updrs2 = df.pheno$updrs2_score,
                            # updrs3 = df.pheno$updrs3_score,
                            # MoCA = df.pheno$moca,
                            col = col_markers, na_col = "grey", 
                            annotation_legend_param = list(legend_direction = "horizontal"))

# ht_meta <- Heatmap(t(trunc.metabolomics[, signatures[[1]]]),
#                   name = "Metabolomics",
#                   row_title = "Metabolomics",
#                   col = bluered(256),
#                   show_column_names = F,
#                   show_row_dend = F,
#                   show_row_names = T,
#                   cluster_columns = F,
#                   # column_split = 3,
#                   column_split = df.pheno$cluster,
#                   column_title = NULL,
#                   top_annotation = col_ha,
#                   height = unit(0.5, "cm"),
#                   heatmap_legend_param = list(legend_direction = "horizontal"))
ht_csf <- Heatmap(t(trunc.csf_proteomics[, signatures[[1]]]),
                  name = "CSF Proteomics",
                  row_title = "CSF Proteomics",
                  col = bluered(256),
                  show_column_names = F,
                  show_row_names = T,
                  show_row_dend = F,
                  cluster_columns = F,
                  cluster_rows = T,
                  column_split = df.pheno$cluster,
                  column_title = NULL,
                  top_annotation = col_ha,
                  height = unit(4, "cm"),
                  row_names_gp = gpar(fontsize = 5),
                  heatmap_legend_param = list(legend_direction = "horizontal"))
ht_urine <- Heatmap(t(trunc.urine_proteomics[, signatures[[2]]]),
                    name = "Urine Proteomics",
                    row_title = "Urine Proteomics",
                    col = bluered(256),
                    show_column_names = F,
                    show_row_names = F,
                    show_row_dend = F,
                    cluster_columns = F,
                    cluster_rows = T,
                    height = unit(5, "cm"),
                    heatmap_legend_param = list(legend_direction = "horizontal"))
ht_rna <- Heatmap(t(trunc.transcriptomics[, signatures[[3]]]),
                  name = "WB Transcriptomics",
                  row_title = "WB Transcriptomics",
                  col = bluered(256),
                  show_column_names = F,
                  show_row_names = F,
                  show_row_dend = F,
                  cluster_columns = F,
                  cluster_rows = T,
                  # Fix heatmap cell size
                  height = unit(4, "cm"),
                  row_names_gp = gpar(fontsize = 6),
                  heatmap_legend_param = list(legend_direction = "horizontal"))

# combine multiple heatmaps
ht_list <- ht_csf %v% ht_urine %v% ht_rna
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/heatmap_gba1_3omics.tiff", res = 300, height = 2000, width = 2000)
draw(ht_list, merge_legends = T, heatmap_legend_side = "right",
     # column_title = "Multi-Omic approach to Neurodegenerative diseases",
     column_title_gp = gpar(fontsize = 20))
dev.off()



# compare clinical data between the two clusters
chisq.test(table(df.pheno$cluster, df.pheno$SEX)) #p-value = 0.5007
chisq.test(table(df.pheno$cluster, df.pheno$PD)) #p-value = 0.3344
chisq.test(table(df.pheno$cluster, df.pheno$CSFSAA_new)) #p-value = 0.5013
t.test(age_at_visit ~ cluster, data = df.pheno) #p-value = 0.05791
t.test(duration_yrs ~ cluster, data = df.pheno) #p-value = 0.2561
t.test(mean_caudate ~ cluster, data = df.pheno) #p-value = 0.3312
t.test(mean_putamen ~ cluster, data = df.pheno) #p-value = 0.471
t.test(mean_striatum ~ cluster, data = df.pheno) #p-value = 0.7237
t.test(log2_ttau ~ cluster, data = df.pheno) #p-value = 0.2426
t.test(log2_ptau ~ cluster, data = df.pheno) #p-value = 0.156
t.test(log2_asyn ~ cluster, data = df.pheno) #p-value = 0.8079
t.test(log2_abeta ~ cluster, data = df.pheno) #p-value = 0.03635
wilcox.test(scopa_gi ~ cluster, data = df.pheno) #p-value = 0.3863
wilcox.test(updrs1_score ~ cluster, data = df.pheno) #p-value = 0.2697
wilcox.test(updrs2_score ~ cluster, data = df.pheno) #p-value = 0.8689
wilcox.test(updrs3_score ~ cluster, data = df.pheno) #p-value = 0.1904
wilcox.test(moca ~ cluster, data = df.pheno) #p-value = 0.9071


# Pathways
library(tidyverse)
library(ggrepel)
library(enrichR) 

setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1")

urine_signatures <- unique(gsub("\\_.*|\\;.*", "", signatures[[2]]))

# Pathways
pathway_enrichr <- function(list, dbs, plot.out, plot.title, table.out){
  res <- enrichr(list, dbs)
  
  # visualize enriched pathways
  tiff(plot.out, res = 300, height = 1500, width = 1500)
  p <- plotEnrich(res[[1]], showTerms = 10, numChar = 50, 
                  y = "Count", orderBy = "P.value",
                  title = plot.title)
  print(p)
  dev.off()
  
  # output tables
  write.table(res[[1]], table.out, sep = '\t', quote = F, row.names = F)
}

# dbs <- listEnrichrDbs()
pathway_enrichr(list = urine_signatures,
                dbs = "KEGG_2021_Human",
                plot.out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/pathways/kegg_urine.tiff",
                plot.title = "GBA1: Urine signatures",
                table.out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/gba1/pathways/kegg_urine.tsv")






library(ggplot2)
temp <- df.pheno %>%
  group_by(cluster, type_mutation) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(percentage = round(n/sum(n),2)) %>%
  as.data.frame()
n <- table(df.pheno$cluster)

tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/proportions_type_mutation.tiff", res = 300, width = 1500, height = 2000)
ggplot(temp, aes(x = cluster, y = percentage, fill = type_mutation))+
  geom_bar(position = "fill", stat = "identity")+
  scale_y_continuous(labels = scales::percent)+
  ylab("Proportion")+
  scale_x_discrete(labels=c("1" = paste0("C1 \n(n = ", n[1], ")"),
                            "2" = paste0("C2 \n(n = ", n[2], ")")))+
  geom_text(aes(label = paste0(percentage*100, "%")),
            position = position_fill(vjust = 0.5))+
  theme_minimal()+
  theme(legend.position = 'none',
                        axis.title = element_text(size = 20),
                        axis.text = element_text(size = 20))
dev.off()

# library(RColorBrewer)
# brewer.pal(n = 8, name = "Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/violin/duration.tiff", res = 300, height = 1500, width = 1500)
ggviolin(df.pheno, x="cluster", y="duration_yrs", color = "cluster",
         add = "boxplot")+
  scale_color_manual(values = c("1"="#D95F02", "2"="#7570B3"))+
  # scale_fill_manual(values = c("1"="#D95F02", "2"="#7570B3"))+
  scale_x_discrete(labels=c("1" = paste0("C1 \n (n = ", n[1], ")"), 
                            "2" = paste0("C2 \n (n = ", n[2], ")")))+
  stat_compare_means(method = "t.test", label.y = 10, size = 6)+
  theme(legend.position = 'none',
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
dev.off()

meta <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_BL.tsv", sep = '\t') %>%
  filter(COHORT == "Healthy Control") %>%
  mutate(cluster = 0) %>%
  select(PATNO, NP2PTOT, cluster)
temp <- df.pheno %>%
  select(PATNO, NP2PTOT, cluster)
df <- rbind(meta, temp)
n <- table(df$cluster)

source("/Volumes/Extreme SSD/work_bidmc/scripts/visualizations/violin.R")
my_comparisons <- list(c("0", "1"), c("0", "2"), c("1", "2"))
violin_func(df = df,
            out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/violin/UPDRS_II_2.tiff", 
            var.x = "cluster", 
            var.y = "NP2PTOT",
            label_y = "UPDRS-II", 
            my_comparisons = my_comparisons)






df.pheno <- merge(df.pheno, urine[, c('PATNO', 'PC1', 'PC2', 'PC3', 'PC4')], by = "PATNO", all.x = T)
df.pheno$cluster <- factor(df.pheno$cluster)
mod <- model.matrix(~ as.factor(cluster)-1+as.factor(SEX)+AGE_AT_VISIT+PC1+PC2, data = df.pheno)
colnames(mod)[1:3] <- c("C1", "C2", "Male")
fit <- lmFit(t(cv.urine_proteomics_top2000[, signatures[[3]]]), mod)
contrast.matrix <- makeContrasts(contrasts = c("C2-C1"), levels = mod)
fitContrasts = contrasts.fit(fit, contrast.matrix)
fit.ebayes <- eBayes(fitContrasts)
n.list <- table(topTable(fit.ebayes, adjust = "BH", n = Inf)$adj.P.Val < 0.05)['TRUE']
n.list
results <- topTable(fit.ebayes, adjust = "BH", n = Inf)
results <- data.frame(symbol = gsub("\\_.*|\\;.*", "", rownames(results)), results)
write.table(results, '/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD_4omics/DE_C2_C1_urine.tsv', sep = '\t', quote = F, row.names = F)

###### box plots of urine protein levels #####

sig.proteins <- data.frame(Cluster = best.cluster.Bayes, REG1A = cv.urine_proteomics_top2000[, "REG1A;REG1B_P05451;P48304"], check.names = F)
n1 <- table(sig.proteins$Cluster)[1]
n2 <- table(sig.proteins$Cluster)[2]

tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD_4omics/urine_REG1A.tiff", res = 300, height = 1500, width = 1500)
ggviolin(sig.proteins, x="Cluster", y="REG1A",
         color="Cluster",
         add = "boxplot")+
  scale_x_discrete(labels=c("1" = paste0("C1 (n = ", n1, ")"), 
                            "2" = paste0("C2 (n = ", n2, ")")))+
  stat_compare_means()
dev.off()


### evaluate the demographic and clinical features
table(df.pheno$cluster, df.pheno$SEX)
table(df.pheno$cluster, df.pheno$type_mutation)
table(df.pheno$cluster, df.pheno$SAA_Status)
table(df.pheno$cluster, df.pheno$APOE)
df.pheno$apoe4 <- ifelse(df.pheno$APOE %in% c("E2/E4", "E3/E4", "E4/E4"), 1, 0)
table(df.pheno$cluster, df.pheno$apoe4)

df.pheno %>% group_by(cluster) %>% 
  summarise(mean(AGE_AT_VISIT, na.rm = T), sd(AGE_AT_VISIT, na.rm = T))

quantile(df.pheno$UPDRS_II[df.pheno$cluster == "1"], probs = c(0.25, 0.5, 0.75), na.rm = T)
quantile(df.pheno$UPDRS_II[df.pheno$cluster == "2"], probs = c(0.25, 0.5, 0.75), na.rm = T)

t.test(df.pheno$AGE_AT_VISIT[df.pheno$cluster == "1"], df.pheno$AGE_AT_VISIT[df.pheno$cluster == "2"])
chisq.test(table(df.pheno$cluster, df.pheno$SAA_Status))

df <- meta %>%
  filter(PATNO %in% comm_id | COHORT == "Healthy Control") %>%
  mutate(cluster = "0",
         UPDRS_I = updrs1_score,
         UPDRS_II = NP2PTOT,
         UPDRS_III = NP3TOT,
         MoCA = MCATOT)
df$cluster[df$PATNO %in% comm_id] <- df.pheno$cluster
df2 <- merge(df, saa, by = "PATNO", all.x = T)

source("/Volumes/Extreme SSD/work_bidmc/scripts/visualizations/violin.R")
setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/violin/2omics")
features <- c("log2_ttau", "log2_ptau", "log2_asyn", "log2_abeta",
             "CAUDATE", "PUTAMEN", "STRIATUM",
             "UPDRS_I", "UPDRS_II", "UPDRS_III", "MoCA",
             "GI")

my_comparisons <- list(c("0", "1"), c("0", "2"), c("1", "2"))
for(f in features){
  violin_func(df = df2, out = paste0(f, ".tiff"),
              var.x = "cluster",
              var.y = f,
              my_comparisons = my_comparisons)
}

t.test(log2_abeta ~ cluster, data = df.pheno)

df.pheno$SAA_Status <- ifelse(df.pheno$SAA_Status == "Positive", 1,  0)
mod <- glm(SAA_Status ~ cluster, data = df.pheno, family = "binomial")
summary(mod)

ord.features <- c("UPDRS_I", "UPDRS_II", "UPDRS_III", "MoCA",
                  "GI")
df_PD <- df.pheno %>%
  mutate(UPDRS_I = factor(UPDRS_I, ordered = T),
         UPDRS_II = factor(UPDRS_II, ordered = T),
         UPDRS_III = factor(UPDRS_III, ordered = T),
         MoCA = factor(MoCA, ordered = T),
         GI = factor(GI, ordered = T))

library(MASS)
ord.test <- function(outcome){
  m <- polr(reformulate("cluster", response = outcome), data = df_PD, Hess = T)
  
  ctable <- coef(summary(m))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  
  ## combined table
  ctable <- cbind(ctable, "p value" = p)
  return(ctable["cluster2", "p value"])
}
for(i in ord.features){
  print(paste0("Outcome: ", i, "; P-value = ", ord.test(outcome = i)))
} 

## GI progression
gi.longit <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_longitudinal.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>% arrange(PATNO, EVENT_ID) %>%
  filter(PATNO %in% comm_id & (EVENT_ID == "BL" | str_detect(EVENT_ID, "^V"))) %>%
  select(PATNO, EVENT_ID, time, GI, NP2PTOT) %>%
  filter(time %in% c(0, 6, 12, 24, 36, 48, 60, 72, 84, 96))

df.longit <- merge(df.pheno[, c("PATNO", "cluster", "duration_yrs", "SEX", "type_mutation")], gi.longit, by = "PATNO")

visit <- factor(gi.longit$EVENT_ID)

tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/longitudinal/longit_GI.tiff")
ggplot(data = df.longit, aes(x = time, y = GI, color = cluster))+
  geom_point(alpha = 0.5)+
  geom_smooth()
  # geom_line(aes(group = cluster))
dev.off()








