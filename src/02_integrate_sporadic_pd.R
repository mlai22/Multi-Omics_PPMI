library(data.table)
library(tidyverse)
library(iClusterPlus)
library(ggpubr)
# library(GenomicRanges)
library(gplots)
# library(lattice) # R graphics

meta <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LONI_phenotype_BL.tsv", sep = '\t')
meta_sPD <- meta %>%
  filter(COHORT == "Parkinson's Disease" & type_mutation == "nomut") # 687 individuals

# CSF Proteomics
setwd('/Volumes/Extreme SSD/work_bidmc/PPMI/proteomics/p151/Proteomic_Analysis')
anno <- read.csv('PPMI_Project_151_pqtl_Analysis_-Table 1.csv')
proteomics_df <- fread('proteomics_220621.tsv', sep = '\t') %>%
  column_to_rownames(var = "PATNO") %>%
  as.data.frame()
probe_human <- anno %>%
  filter(ORGANISM == "Human") %>%
  select(SOMA_SEQ_ID)
proteomics_human <- proteomics_df[, colnames(proteomics_df) %in% probe_human$SOMA_SEQ_ID] #1158x4779
name_symbol <- anno[match(colnames(proteomics_human), anno[, 'SOMA_SEQ_ID']), 'TARGET_GENE_SYMBOL'] #4779 SOMAmers targeting 4102 unique proteins
colnames(proteomics_human) <- name_symbol

# Urine Proteomics
urine <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/Urine_Proteomics/DF2_log2_pc4.tsv", sep = '\t') %>%
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
setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/metabolomics_PPMI180/Metabolomic_Analysis")
dat5 <- fread("Metabolomic_Analysis_of_LRRK2_PD_5_of_5.csv") %>%
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
rna <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/PPMI_rnaseq/log2CPM_top2000_chr22_coding.tsv", sep = '\t')
rna2 <- t(rna[, -1])
colnames(rna2) <- rna$Gene
rownames(rna2) <- gsub("PP-", "", rownames(rna2))
rna2 <- rna2[order(as.numeric(rownames(rna2))),]

id <- rownames(proteomics_human)[rownames(proteomics_human) %in% pheno4_sPD$PATNO] #377
# id <- rownames(urine2_qc)[rownames(urine2_qc) %in% id] #228
id <- rownames(rna2)[rownames(rna2) %in% id] #212
id <- rownames(metabolics_qc)[rownames(metabolics_qc) %in% id] #69

df.metabolics_sPD <- as.matrix(metabolics_qc[rownames(metabolics_qc) %in% id,]) #69x281
df.proteomics_sPD <- as.matrix(proteomics_human[rownames(proteomics_human) %in% id,]) #69x4779
# df.urine_proteomics_sPD <- as.matrix(urine2_qc[rownames(urine2_qc) %in% id,]) #69x2170
df.transcriptomics_sPD <- as.matrix(rna2[rownames(rna2) %in% id,]) #69x2000

# most variable features
compute_cv <- function(x){ sd(x)/mean(x) }
cv.proteomics_sPD <- apply(df.proteomics_sPD, 2, compute_cv)
cv.proteomics_sPD <- sort(cv.proteomics_sPD, decreasing = T)
cv.proteomics_sPD_top2000 <- names(head(cv.proteomics_sPD, 2000))
cv.proteomics_sPD_top2000 <- df.proteomics_sPD[, cv.proteomics_sPD_top2000]
colnames(cv.proteomics_sPD_top2000) <- make.unique(colnames(cv.proteomics_sPD_top2000), sep = '.')

# cv.urine_proteomics_sPD <- apply(df.urine_proteomics_sPD, 2, compute_cv)
# cv.urine_proteomics_sPD <- sort(cv.urine_proteomics_sPD, decreasing = T)
# cv.urine_proteomics_sPD_top2000 <- names(head(cv.urine_proteomics_sPD, 2000))
# cv.urine_proteomics_sPD_top2000 <- df.urine_proteomics_sPD[, cv.urine_proteomics_sPD_top2000]

################################################################################
#                                iClusterBayes                                 #
################################################################################
set.seed(123)
date()
bayfit <- tune.iClusterBayes(cpus = 12, 
                             dt1 = df.metabolics_sPD, dt2 = cv.proteomics_sPD_top2000, dt3 = df.transcriptomics_sPD,
                             type = c('gaussian', 'gaussian', 'gaussian'),
                             K = 1:6, n.burnin = 18000, # number of MCMC burnin
                             n.draw = 22000, # number of MCMC draw
                             prior.gamma = c(0.1, 0.1, 0.1),
                             sdev = 0.05, pp.cutoff = 0.5)
date()
save.image(file = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_meta_csf_rna_sPD/PPMI_meta_csf_rna_sPD_coding2.RData")
load("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_4omics_sPD/PPMI_4omics_sPD_coding.RData")

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
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_4omics_sPD/BIC_deviance_plot_4omics_sPD2.tiff", res = 300, height = 1500, width = 1500)
par(mar=c(4.0, 4.0, 0.5, 0.5), mfrow=c(1, 2))
plot(1:nK, allBIC, type="b", xlab="k", ylab="BIC", pch=c(19, 1, 1, 1, 1, 1))
plot(1:nK, devratio, type="b", xlab="k", ylab="Deviance ratio", pch=c(19, 1, 1, 1, 1, 1))
dev.off()

k <- 1
##### compare the results by iClusterPlus and iClusterBayes
best.cluster.Bayes = bayfit$fit[[k]]$clusters
table(best.cluster.Bayes)
# best.cluster is the clusters generated by iClusterPlus; See the above code
# table(best.cluster, best.cluster.Bayes)

df.pheno <- pheno_sPD[pheno_sPD$PATNO %in% rownames(cv.proteomics_sPD_top2000),]
df.pheno <- df.pheno[order(as.numeric(df.pheno$PATNO)),]
rownames(df.pheno) <- df.pheno$PATNO
df.pheno$cluster <- as.character(best.cluster.Bayes)

df_HC <- data.frame(pheno[pheno$COHORT == "Control",], cluster = "0") #186
df_HC_sPD <- rbind(df_HC, df.pheno)
df_HC_sPD2 <- merge(df_HC_sPD, pheno_BL[, c("PATNO", "log2_ttau", "log2_ptau", "log2_asyn", "log2_abeta",
                                            "updrs1_score", "NP2PTOT", "NP3TOT", "MCATOT",
                                            "CAUDATE", "PUTAMEN", "STRIATUM")],
                    by = "PATNO")
################################################################################
#                               Generate Heatmap                               #
################################################################################
# # Ggenomic features with posterior probability > 0.5 were selected
# col.scheme = alist()
# col.scheme[[1]] = bluered(256)
# col.scheme[[2]] = bluered(256)
# col.scheme[[3]] = bluered(256)
# col.scheme[[4]] = bluered(256)

# scale each dataset
scale <- function(x){ (x-mean(x))/sd(x) }
scale.metabolics_sPD <- apply(df.metabolics_sPD, 2, scale)
scale.metabolics_sPD[scale.metabolics_sPD < -4] <- -4
scale.metabolics_sPD[scale.metabolics_sPD > 4] <- 4
scale.proteomics_sPD <- apply(cv.proteomics_sPD_top2000, 2, scale)
scale.proteomics_sPD[scale.proteomics_sPD < -2] <- -2
scale.proteomics_sPD[scale.proteomics_sPD > 2] <- 2
# scale.urine_proteomics_sPD <- apply(cv.urine_proteomics_sPD_top2000, 2, scale)
# scale.urine_proteomics_sPD[scale.urine_proteomics_sPD < -5] <- -5
# scale.urine_proteomics_sPD[scale.urine_proteomics_sPD > 5] <-5
scale.transcriptomics_sPD <- apply(df.transcriptomics_sPD, 2, scale)
scale.transcriptomics_sPD[scale.transcriptomics_sPD < -5] <- -5
scale.transcriptomics_sPD[scale.transcriptomics_sPD > 5] <- 5

# tiff("/Users/mei-yulai/Desktop/iCluster_PPMI/PPMI_meta_csf_rna_sPD.tiff", res = 300, height = 2000, width = 2000)
# plotHMBayes(fit=bayfit$fit[[k]], 
#             datasets=list(scale.metabolics_sPD, scale.proteomics_sPD, scale.transcriptomics_sPD),
#             type=c("gaussian", "gaussian", "gaussian"), col.scheme = col.scheme,
#             threshold=c(0.5, 0.5, 0.5), row.order=c(T, T, T), scale=c("none", "none", "none"), 
#             plot.chr=NULL, sparse=c(T, T, T), cap=NULL)
# dev.off()

features <- alist()
features[[1]] <- colnames(df.metabolics_sPD)
features[[2]] <- colnames(cv.proteomics_sPD_top2000)
# features[[3]] <- colnames(cv.urine_proteomics_sPD_top2000)
features[[3]] <- colnames(df.transcriptomics_sPD)

signatures <- alist()
for(i in 1:3){
  signatures[[i]] <- (features[[i]])[which(bayfit$fit[[k]]$beta.pp[[i]] > 0.5)]
}
length(signatures[[1]])
length(signatures[[2]])
length(signatures[[3]])
# length(signatures[[4]])

### venn diagram ###
library(VennDiagram)
csf <- signatures[[2]]
# urine <- signatures[[3]]
# urine <- gsub("\\_.*|\\;.*", "", urine)
rna <- signatures[[3]]

v <- venn.diagram(x = list(csf, rna),
                  category.names = c('CSF_proteomics', 'WB_rna'),
                  filename = NULL,
                  cex=1.5,
                  cat.pos = c(135, 180),
                  cat.dist = c(0.01, 0.025),
                  margin = 0.05,
                  col = c("#440154ff", "#F39200FF"),
                  fill = c("#440154ff", "#F39200FF"))
# col = c("#440154ff", "#31B7BC99", "#F39200FF"),
# fill = c("#440154ff", "#31B7BC99", "#F39200FF"))
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_meta_csf_rna_sPD/venn_coding_meta_csf_rna.tiff", res = 300, width = 1500, height = 1500)
grid.newpage()
grid.draw(v)
dev.off()


# Heatmaps
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

ttau <- (df.pheno$log2_ttau - mean(df.pheno$log2_ttau, na.rm = T))/sd(df.pheno$log2_ttau, na.rm = T)
ptau <- (df.pheno$log2_ptau - mean(df.pheno$log2_ptau, na.rm = T))/sd(df.pheno$log2_ptau, na.rm = T)
asyn <- (df.pheno$log2_asyn - mean(df.pheno$log2_asyn, na.rm = T))/sd(df.pheno$log2_asyn, na.rm = T)
abeta <- (df.pheno$log2_abeta - mean(df.pheno$log2_abeta, na.rm = T))/sd(df.pheno$log2_abeta, na.rm = T)


col_markers <- list(subtype = c("1" = "#7570B3", "2" = "#1B9E77", "3" = "#F39200FF"),
                    ttau = colorRamp2(range(ttau, na.rm = T), colors = c("white", "darkgreen")),
                    ptau = colorRamp2(range(ptau, na.rm = T), colors = c("white", "darkgreen")),
                    asyn = colorRamp2(range(asyn, na.rm = T), colors = c("white", "darkgreen")),
                    abeta = colorRamp2(range(abeta, na.rm = T), colors = c("white", "darkgreen")))
col_ha <- HeatmapAnnotation(subtype = df.pheno$cluster, 
                            ttau = ttau, ptau = ptau, asyn = asyn, abeta = abeta,
                            col = col_markers, na_col = "black", 
                            annotation_legend_param = list(legend_direction = "horizontal"))

ht_meta <- Heatmap(t(scale.metabolics_sPD[, signatures[[1]]]),
                   name = "Metabolomics",
                   row_title = "Metabolomics",
                   col = bluered(256),
                   show_column_names = F,
                   show_row_dend = F,
                   show_row_names = T,
                   cluster_columns = T,
                   # column_split = 3,
                   column_split = df.pheno$cluster,
                   column_title = NULL,
                   top_annotation = col_ha,
                   height = unit(6, "cm"),
                   heatmap_legend_param = list(legend_direction = "horizontal"))
ht_csf <- Heatmap(t(scale.proteomics_sPD[, signatures[[2]]]),
                  name = "CSF Proteomics",
                  row_title = "CSF Proteomics",
                  col = bluered(256),
                  show_column_names = F,
                  show_row_dend = F,
                  cluster_columns = F,
                  column_split = df.pheno$cluster,
                  height = unit(6, "cm"),
                  row_names_gp = gpar(fontsize = 6),
                  heatmap_legend_param = list(legend_direction = "horizontal"))
# ht_urine <- Heatmap(t(scale.urine_proteomics_sPD[, signatures[[3]]]),
#                     name = "Urine Proteomics",
#                     row_title = "Urine Proteomics",
#                     col = bluered(256),
#                     show_column_names = F,
#                     show_row_names = F,
#                     show_row_dend = F,
#                     cluster_columns = F,
#                     height = unit(4, "cm"),
#                     heatmap_legend_param = list(legend_direction = "horizontal"))
ht_rna <- Heatmap(t(scale.transcriptomics_sPD[, signatures[[3]]]),
                  name = "WB Transcriptomics",
                  row_title = "WB Transcriptomics",
                  col = bluered(256),
                  show_column_names = F,
                  show_row_names = F,
                  show_row_dend = F,
                  cluster_columns = F,
                  # Fix heatmap cell size
                  height = unit(6, "cm"),
                  row_names_gp = gpar(fontsize = 6),
                  heatmap_legend_param = list(legend_direction = "horizontal"))

# combine multiple heatmaps
ht_list <- ht_meta %v% ht_csf %v% ht_rna
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_meta_csf_rna_sPD/heatmap_meta_csf_rna_coding2.tiff", res = 300, height = 3000, width = 3000)
draw(ht_list, merge_legends = T, heatmap_legend_side = "right")
dev.off()


# ### features with most difference
# diff <- NULL
# temp <- df.transcriptomics_PD[, which(bayfit$fit[[k]]$beta.pp[[3]] > 0.5)]
# for(i in 1:ncol(temp)){
#   m1 <- mean(temp[best.cluster.Bayes == 1, i])
#   m2 <- mean(temp[best.cluster.Bayes == 2, i])
#   diff[i] <- abs(m1 - m2)
# }
# names(diff) <- colnames(temp)
# head(sort(diff, decreasing = T), 10)
# 
# 
# ### row-ordered features
# image.data <- df.transcriptomics_PD[, which(bayfit$fit[[k]]$beta.pp[[3]] > 0.5)]
# diss = 1 - cor(image.data, use = "na.or.complete")
# hclust.fit = hclust(as.dist(diss))
# gorder = hclust.fit$order
# image.data = image.data[, gorder]
# 

################################################################################
#                                   DE analysis                                #
################################################################################
# library(limma)
# 
# df.pheno <- merge(df.pheno, urine[, c('PATNO', 'PC1', 'PC2', 'PC3', 'PC4')], by = "PATNO", all.x = T)
# df.pheno$cluster <- factor(df.pheno$cluster)
# mod <- model.matrix(~ as.factor(cluster)-1+as.factor(SEX)+AGE_AT_VISIT+PC1+PC2, data = df.pheno)
# colnames(mod)[1:3] <- c("C1", "C2", "Male")
# fit <- lmFit(t(cv.urine_proteomics_sPD_top2000[, signatures[[3]]]), mod)
# contrast.matrix <- makeContrasts(contrasts = c("C2-C1"), levels = mod)
# fitContrasts = contrasts.fit(fit, contrast.matrix)
# fit.ebayes <- eBayes(fitContrasts)
# n.list <- table(topTable(fit.ebayes, adjust = "BH", n = Inf)$adj.P.Val < 0.05)['TRUE']
# n.list
# results <- topTable(fit.ebayes, adjust = "BH", n = Inf)
# results <- data.frame(symbol = gsub("\\_.*|\\;.*", "", rownames(results)), results)
# write.table(results, '/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_4omics_sPD/DE_C2_C1_urine.tsv', sep = '\t', quote = F, row.names = F)

# model for transcriptomics
library(DESeq2)
df_HC <- data.frame(pheno4[pheno4$COHORT == "Control" & !is.na(pheno4$COHORT),], cluster = "0")
df_HC_sPD <- rbind(df.pheno, df_HC)

rna_count <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/PPMI_rnaseq/PPMI_BL_rna_chr22_coding.tsv", sep = '\t', header = T)
rna2 <- as.data.frame(rna_count[, -c(1:3)])
rownames(rna2) <- make.unique(rna_count$hgnc_symbol, sep = ".")
colnames(rna2) <- gsub("-BLM0T1", "", colnames(rna2))
colnames(rna2) <- gsub("PP-", "", colnames(rna2))
rna2 <- rna2[, order(as.numeric(names(rna2)))]
rna2 <- rna2[rownames(rna2) %in% signatures[[3]], colnames(rna2) %in% df_HC_sPD$PATNO]

df_HC_sPD <- df_HC_sPD[df_HC_sPD$PATNO %in% colnames(rna2),]
df_HC_sPD <- df_HC_sPD[order(as.numeric(df_HC_sPD$PATNO)),]
rownames(df_HC_sPD) <- df_HC_sPD$PATNO
all(colnames(rna2) == rownames(df_HC_sPD))
df_HC_sPD$cluster <- factor(df_HC_sPD$cluster, levels = c("0", "1", "2", "3"))

dds <- DESeqDataSetFromMatrix(countData = rna2,
                              colData = df_HC_sPD,
                              design = ~ cluster+age_at_baseline+sex)
# pre-filtering: removing rows in which there are very few reads
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]

library("BiocParallel")
register(MulticoreParam(4))

de_dds <- DESeq(dds)
res <- results(de_dds, contrast=c("cluster", "1", "2"), 
               alpha = 0.05, lfcThreshold = 0.1)
res <- res[order(res$pvalue),]
summary(res)
write.table(cbind(Gene=rownames(res), res), "/Users/meiyu/Desktop/work_bidmc/PPMI/iCluster/DE_top25/rna_C1_HC_6.tsv",
            sep = '\t', quote = F, row.names = F)


###### box plots of urine protein levels #####
sig.proteins <- data.frame(Cluster = best.cluster.Bayes, PSMA7 = cv.urine_proteomics_sPD_top2000[, "PSMA7_O14818;O14818-2"], check.names = F)
HC <- data.frame(Cluster = 0, PSMA7 = urine[urine$PATNO %in% pheno4$PATNO[pheno4$COHORT == "Control"], "PSMA7_O14818;O14818-2"], check.names = F)
colnames(HC) <- colnames(sig.proteins)
sig.proteins <- rbind(sig.proteins, HC)
n0 <- table(sig.proteins$Cluster)[1]
n1 <- table(sig.proteins$Cluster)[2]
n2 <- table(sig.proteins$Cluster)[3]

my_comparisons <- list(c("0", "1"), c("0", "2"), c("1", "2"))
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_4omics_sPD/urine_PSMA7.tiff", res = 300, height = 1500, width = 1500)
ggviolin(sig.proteins, x="Cluster", y="PSMA7",
         color="Cluster",
         add = "boxplot")+
  scale_x_discrete(labels=c("0" = paste0("HC (n = ", n0, ")"),
                            "1" = paste0("C1 (n = ", n1, ")"), 
                            "2" = paste0("C2 (n = ", n2, ")")))+
  stat_compare_means(comparisons = my_comparisons)
dev.off()


### evaluate the demographic and clinical features
table(df_HC_sPD$cluster, df_HC_sPD$sex)
table(df_HC_sPD$cluster, df_HC_sPD$APOE_genotype)
df_HC_sPD %>% group_by(cluster) %>% 
  summarise(mean(log2_abeta, na.rm = T), sd(log2_abeta, na.rm = T))

df <- df_HC_sPD2[, c('MCATOT', 'cluster')]
df$cluster <- factor(df$cluster, levels = c("0", "1", "2"))
df <- df[complete.cases(df),]
n <- table(df$cluster)

my_comparisons <- list(c("0", "1"), c("0", "2"), c("1", "2"))
tiff("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_4omics_sPD/MCATOT.tiff", res = 300, height = 1500, width = 1500)
ggviolin(df, x="cluster", y="MCATOT", 
         color="cluster",
         add = "boxplot")+
  scale_x_discrete(labels=c("0" = paste0("HC \n (n = ", n[1], ")"),
                            "1" = paste0("C1 \n (n = ", n[2], ")"),
                            "2" = paste0("C2 \n (n = ", n[3], ")")))+
  stat_compare_means(comparisons = my_comparisons)
# stat_compare_means(method = "t.test", label.y = 5000)
dev.off()

