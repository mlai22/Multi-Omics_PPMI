library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice) # R graphics
library(data.table)
library(tidyverse)

setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/coding_top2000_meta_csf_rna_sPD")
load("PPMI_meta_csf_rna_sPD_coding.RData")

k <- 2
##### compare the results by iClusterPlus and iClusterBayes
best.cluster.Bayes = bayfit$fit[[k]]$clusters
table(best.cluster.Bayes)

##### top features #####
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





##### phenotype #####
pheno <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/pheno/Bruno/PPMI_New_Pheno_080122.csv")
pheno_sPD <- pheno[(pheno$COHORT == 'PD' & pheno$GBA_PATHVAR == 0 & pheno$PATNO %in% rownames(cv.proteomics_sPD_top2000)) | pheno$COHORT == "Control",] #255 

pheno2 <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/pheno/PPMI_pheno_BL_amp.tsv", sep = '\t')
pheno2$PATNO <- gsub("PP-", "", pheno2$participant_id)
asyn <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/pheno/Bruno/PPMI_2022-0001_Serrano-Fernandez_snapshot.csv")
pheno2 <- merge(pheno2, asyn[, c('PATNO', 'asyn')], by = 'PATNO', all.x = T)
pheno2$mean_striatum <- (pheno2$mean_caudate + pheno2$mean_putamen)/2
pheno2_sPD <- pheno2[pheno2$PATNO %in% pheno_sPD$PATNO,] #288

df.pheno <- merge(pheno_sPD, pheno2_sPD, by = "PATNO")
df.pheno$cluster <- "0"
df.pheno$cluster[df.pheno$PATNO %in% rownames(cv.proteomics_sPD_top2000)] <- as.character(best.cluster.Bayes)
df.pheno$cluster <- recode_factor(df.pheno$cluster, `0` = "HC", `1` = "C1", '2' = "C2", '3' = "C3")

df.pheno$log2_ttau <- log2(df.pheno$Tau+1)
df.pheno$log2_ptau <- log2(df.pheno$p.Tau+1)
df.pheno$log2_asyn <- log2(df.pheno$asyn+1)
df.pheno$log2_abeta <- log2(df.pheno$Abeta+1)

##### Demographics #####
# sex
table(df.pheno$sex, df.pheno$cluster)
#age_at_BL
df.pheno %>% group_by(cluster) %>% summarise(mean(age_at_baseline), sd(age_at_baseline))
#mutations
table(df.pheno$type_mutation, df.pheno$cluster)
table(df.pheno$APOE_genotype, df.pheno$cluster)


# merge phenotype with proteomic PCs
proteomics_PC <- read.table("/Users/meiyu/Desktop/work_bidmc/PPMI/pheno/Bruno/Proteomics_PC.tsv", sep = '\t', header = T)
pheno_sub[pheno_sub$PATNO %in% proteomics_PC$PATNO, c('PC1', 'PC2', 'PC3', 'PC4')] <- proteomics_PC[proteomics_PC$PATNO %in% pheno_sub$PATNO, -1]
# extract proteomcis, metabolics, and transcriptomics
setwd('/Users/meiyu/Downloads')
anno <- read.csv('PPMI_Project_151_pqtl_Analysis_-Table 1.csv')
proteomics_df <- read.csv('proteomics_220621.tsv', sep = '\t', check.names = F, header = T)
proteomics <- proteomics_df[, -1]
matched_probe <- anno[anno$SOMA_SEQ_ID %in% colnames(proteomics),]
probe_nonhuman <- matched_probe[matched_probe$ORGANISM != "Human", 'SOMA_SEQ_ID']
proteomics_human <- proteomics[, !(colnames(proteomics) %in% probe_nonhuman)] #1158x4779
rownames(proteomics_human) <- proteomics_df[, 1]
name_symbol <- anno[match(colnames(proteomics_human), anno[, 'SOMA_SEQ_ID']), 'TARGET_GENE_SYMBOL'] #4779 SOMAmers targeting 4102 unique proteins
df.proteomics <- proteomics_human[rownames(proteomics_human) %in% pheno_sub$PATNO,]
pheno.proteomics <- pheno_sub[rownames(pheno_sub) %in% rownames(df.proteomics),]
table(pheno.proteomics$cohort)

setwd("/Users/meiyu/Desktop/work_bidmc/PPMI/metabolomics_PPMI180/Metabolomic_Analysis")
dat5 <- fread("Metabolomic_Analysis_of_LRRK2_PD_5_of_5.csv")
dat5 <- dat5[order(dat5$PATNO),]
# select UNITS:adjusted_area_ratio
dat5 <- dat5 %>% filter(UNITS == 'adjusted_area_ratio')
dat5 <- dat5[, c('PATNO', 'SEX', 'COHORT', 'CLINICAL_EVENT', 'TESTNAME', 'TESTVALUE')]
# reshape the data from long to wide
dat5_wide <- spread(dat5, TESTNAME, TESTVALUE)
# unique PATNO ID
dat5_uniq <- dat5_wide[!duplicated(dat5_wide$PATNO),] #607
metabolics <- dat5_uniq[, -c(1:4)] #607x298
rownames(metabolics) <- dat5_uniq$PATNO
# QC: exclude metabolites with > 20% missing values
count <- apply(metabolics, 2, function(x){
  sum(is.na(x))
})
metabolics_qc <- metabolics[, -which(count > nrow(metabolics)*0.2)] #607x281
metabolics_qc <- log2(metabolics_qc)

df.metabolics <- metabolics_qc[rownames(metabolics_qc) %in% rownames(pheno_sub),]
pheno.metabolics <- pheno_sub[rownames(pheno_sub) %in% rownames(metabolics_qc),]
# add caffeine in the model of metabolics
ca <- read.csv("/Users/meiyu/Desktop/amp-pd/2022_v3release_1115/clinical/Caffeine_history.csv")
ca_PP <- ca[grep("PP-", ca$participant_id),]
ca_PP$PATNO <- as.numeric(gsub("PP-", "", ca_PP$participant_id))
ca_PP <- ca_PP[order(ca_PP$PATNO),]
rownames(ca_PP) <- ca_PP$PATNO
pheno.metabolics$caffeine[rownames(pheno.metabolics) %in% rownames(ca_PP)] <- ca_PP$caff_drinks_current_use[rownames(ca_PP) %in% rownames(pheno.metabolics)]
table(pheno.metabolics$cohort)

rna <- fread("/Users/meiyu/Desktop/work_bidmc/PPMI/PPMI_rnaseq/log2CPM_top25.tsv", sep = '\t', header = T)
rna2 <- t(rna[, -1])
colnames(rna2) <- rna$Gene
rownames(rna2) <- gsub("PP.", "", rownames(rna2))
rna2 <- rna2[order(as.numeric(rownames(rna2))),]
df.rna <- rna2[rownames(rna2) %in% rownames(pheno_sub),]
pheno.rna <- pheno_sub[rownames(pheno_sub) %in% rownames(df.rna),]
table(pheno.rna$cohort)

# model for proteomics
mod <- model.matrix(~ as.factor(cluster)-1+as.factor(sex)+age_at_BL+PC1+PC2+PC3+PC4, data = pheno.proteomics)
colnames(mod)[1:4] <- c("HC", "C1", "C2", "Male")
fit <- lmFit(t(df.proteomics), mod)
contrast.matrix <- makeContrasts(contrasts = c("C1-HC"), levels = mod)
fitContrasts = contrasts.fit(fit, contrast.matrix)
fit.ebayes <- eBayes(fitContrasts)
n.list <- table(topTable(fit.ebayes, adjust = "BH", n = Inf)$adj.P.Val < 0.05)['TRUE']
n.list
results <- topTable(fit.ebayes, adjust = "BH", n = Inf)
results <- data.frame(symbol = anno[match(row.names(results), anno[, 'SOMA_SEQ_ID']), 'TARGET_GENE_SYMBOL'],
                      probe=rownames(results), results)
write.table(results, '/Users/meiyu/Desktop/work_bidmc/PPMI/iCluster/DE_top25/proteomics_C2_HC_13.tsv', sep = '\t', quote = F, row.names = F)

# model for metabolics
mod <- model.matrix(~ as.factor(cluster)-1+as.factor(sex)+age_at_BL+Sum_LEDD, data = pheno.metabolics)
colnames(mod)[1:4] <- c("HC", "C1", "C2", "Male")
fit <- lmFit(t(df.metabolics), mod)
contrast.matrix <- makeContrasts(contrasts = c("C2-HC"), levels = mod)
fitContrasts = contrasts.fit(fit, contrast.matrix)
fit.ebayes <- eBayes(fitContrasts)
n.list <- table(topTable(fit.ebayes, adjust = "BH", n = Inf)$adj.P.Val < 0.05)['TRUE']
n.list
results <- topTable(fit.ebayes, adjust = "BH", n = Inf)
results <- data.frame(symbol = rownames(results), results)
write.table(results, '/Users/meiyu/Desktop/work_bidmc/PPMI/iCluster/DE_top25/metabolics_C2_HC_1.tsv', sep = '\t', quote = F, row.names = F)





###### bar plot #####
# APOE4 protportions
library(scales) # convert float to percentage
df <- pheno_sub[, c('PATNO', 'cluster', 'APOE')]
df$cluster <- factor(df$cluster, levels = c('0', "1", "2"))
df$APOE <- factor(df$APOE, levels = c("E2/E2", "E2/E3", "E3/E3", "E2/E4", "E3/E4", "E4/E4", "NA"))
df <- df %>%
  group_by(cluster, .drop = F) %>% 
  dplyr::count(APOE) %>% 
  mutate(percentage = label_percent()( round(n/sum(n), 2) ))

tiff("/Users/meiyu/Desktop/work_bidmc/PPMI/iCluster/DE_top25/apoe_proportions.tiff", res = 300, height = 1800, width = 2000)
ggbarplot(df, x="cluster", y="n",
          fill="APOE", 
          label = df$percentage, lab.size = 3, 
          palette = c("#007B3D99", "#95C11F99", "#31B7BC99", "#F3920099", "#F39200FF", "#D5131799", "#86868699"), 
          position = position_dodge(0.8))+
  scale_x_discrete(labels=c('HC\n(n = 230)', 
                            'C1\n(n = 108)', 
                            'C2\n(n = 112)'))
  # stat_compare_means(comparisons = my_comparisons)+
  # stat_compare_means(label.y = 150)
dev.off()

