library(tidyverse)
library(data.table)
library(DESeq2)

df_pheno <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/df_PD3omics.tsv", sep = '\t') %>%
  mutate(cluster = factor(cluster),
         type_mutation = factor(type_mutation, levels = c('nomut', 'GBA', 'LRRK2'))) %>%
  select(PATNO, cluster, SEX, AGE_AT_VISIT, duration_yrs, type_mutation) %>%
  arrange(PATNO) %>%
  column_to_rownames("PATNO")

################################################################################
#                                   DE analysis                                #
################################################################################
source("/Volumes/Extreme SSD/work_bidmc/scripts/DE/DE_limma.R")
source("/Volumes/Extreme SSD/work_bidmc/scripts/visualizations/volcano.R")
source("/Volumes/Extreme SSD/work_bidmc/scripts/visualizations/volcano_rna.R")

###############################################################################
#                               Transcriptomics                               #
###############################################################################
setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/PPMI_rnaseq")
rna_count <- fread("PPMI_BL_rna_chr22_coding.tsv", sep = '\t')
df_rna <- rna_count %>%
  mutate(hgnc_symbol = make.unique(hgnc_symbol)) %>%
  column_to_rownames("hgnc_symbol") %>%
  select(-GENEID, -Geneid)
sample <- gsub("-BLM0T1", "", colnames(df_rna))
df_rna <- df_rna[, which(sample %in% rownames(df_pheno))]
df_rna <- df_rna[, order(colnames(df_rna))] # 645
colnames(df_rna) <- gsub("-BLM0T1", "", colnames(df_rna))

###################
#     DESeq2      #
###################
dds <- DESeqDataSetFromMatrix(countData = df_rna,
                              colData = df_pheno,
                              design = ~ cluster+SEX+AGE_AT_VISIT+duration_yrs+type_mutation)
# pre-filtering: removing rows in which there are very few reads
keep <- rowSums(counts(dds) >= 5) >= 10
dds <- dds[keep,]

# speed up
library("BiocParallel")
register(MulticoreParam(4))

dds <- DESeq(dds)

res <- results(dds, contrast = c("cluster", "1", "2"), alpha = 0.05)
summary(res)
# res2 <- lfcShrink(dds, coef = "clusters_B_vs_A")
# summary(res2)
results <- res[order(res$padj),]
results <- data.frame(Gene=rownames(results), results)
write.table(results, "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_transcriptomics.tsv", 
            sep = '\t', quote = F, row.names = F)

# volcano plot
results <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_transcriptomics.tsv", sep = '\t')
volcano_rna(de = results, lfc = 0, pval = 0.05, 
            output = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_transcriptomics.tiff",
            title = "WB transcriptomics")

###############################################################################
#                             End: Transcriptomics                            #
###############################################################################

###############################################################################
#                                Metabolomics                                 #
###############################################################################
setwd("/Volumes/Extreme SSD/work_bidmc/PPMI/metabolomics_PPMI180/Metabolomic_Analysis")

dat5 <- fread("Metabolomic_Analysis_of_LRRK2_PD_5_of_5.csv")
dat5 <- as.data.frame(dat5)
# select UNITS:adjusted_area_ratio
metabolomics <- dat5 %>% filter(UNITS == 'adjusted_area_ratio') %>% 
  arrange(PATNO, CLINICAL_EVENT) %>%
  select(PATNO, CLINICAL_EVENT, TESTNAME, TESTVALUE) %>%
  spread(TESTNAME, TESTVALUE) %>% # reshape the data from long to wide
  distinct(PATNO, .keep_all = T) %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>%
  arrange(PATNO) %>%
  column_to_rownames("PATNO") #607

df_metabolomics <- metabolomics[, -1] #607x298
# QC: exclude metabolites with > 20% missing values
count <- apply(df_metabolomics, 2, function(x){
  sum(is.na(x))
})
metabolics_qc <- df_metabolomics[, -which(count > nrow(df_metabolomics)*0.2)] #607x281
metabolics_qc <- log2(metabolics_qc)

# update the age at visit due to the different visit time
df.pheno <- df_pheno[rownames(df_pheno) %in% rownames(metabolics_qc),] #176
metabolics_qc <- metabolics_qc[rownames(metabolics_qc) %in% rownames(df.pheno),]

ledd <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/meta_data/LEDD_BL_final.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO))
df.pheno <- merge(df.pheno, ledd, by.x = "row.names", by.y = "PATNO") %>%
  column_to_rownames("Row.names")

mod <- model.matrix(~ cluster-1 + AGE_AT_VISIT + SEX + duration_yrs + type_mutation + Sum_LEDD, data = df.pheno)
colnames(mod) <- c('C1', 'C2', 'AGE_AT_VISIT', 'Male', 'duration_yrs', 'GBA1', 'LRRK2', 'Sum_LEDD')
res <- de_func(model = mod, feature = metabolics_qc, contrasts = "C1-C2", 
               table.out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_metabolite.tsv")

###############################################################################
#                             End: Metabolomics                               #
###############################################################################

###############################################################################
#                               Urine Proteomics                              #
###############################################################################
urine_proteomics <- fread("/Volumes/Extreme SSD/work_bidmc/PPMI/proteomics/Urine_Proteomics/DF2_log2_pc4.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>% 
  arrange(PATNO) %>%
  as.data.frame()
df_proteomics <- as.data.frame(urine_proteomics[, -c(1:9)])
rownames(df_proteomics) <- urine_proteomics$PATNO

df.pheno <- merge(urine_proteomics[, c(1, 6:9)], df_pheno, by.x = "PATNO", by.y = "row.names") %>%
  column_to_rownames("PATNO") #373
df_proteomics <- df_proteomics[rownames(df_proteomics) %in% rownames(df.pheno),]

mod <- model.matrix(~ cluster-1 + AGE_AT_VISIT + SEX + duration_yrs + type_mutation + PC1 + PC2 + PC3 + PC4, data = df.pheno)
colnames(mod)[1:7] <- c('C1', 'C2', 'AGE_AT_VISIT', 'Male', 'duration_yrs', 'GBA1', 'LRRK2')
res <- de_func(model = mod, feature = df_proteomics, contrasts = "C1-C2", 
               table.out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_proteomics_urine.tsv")

results <- read.csv("/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_proteomics_urine.tsv", sep = '\t')
results$Symbol <- gsub(".*_", "", results$Symbol)
results$Symbol <- gsub("\\;.*", "", results$Symbol)

volcano_func(df = results, thres_fc = 0, pval = results$adj.P.Val, 
             output = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_proteomics_urine2.tiff", 
             title = "Urine proteomics")

###############################################################################
#                             End: Urine Proteomics                           #
###############################################################################

###############################################################################
#                               CSF Proteomics                                #
###############################################################################
setwd('/Volumes/Extreme SSD/work_bidmc/PPMI/proteomics/p151/Proteomic_Analysis')
anno <- read.csv('PPMI_Project_151_pqtl_Analysis_-Table 1.csv')
proteomics <- fread('proteomics_220621.tsv', sep = '\t', check.names = F, header = T) %>%
  as.data.frame() %>%
  mutate(PATNO = paste0("PP-", PATNO)) %>%
  arrange(PATNO) %>%
  column_to_rownames("PATNO")

probe_human <- anno %>% 
  filter(SOMA_SEQ_ID %in% colnames(proteomics) & ORGANISM == "Human") %>% 
  distinct(SOMA_SEQ_ID)
proteomics_human <- proteomics[, probe_human[, 1]] #1158x4779
name_symbol <- anno[match(colnames(proteomics_human), anno[, 'SOMA_SEQ_ID']), 'TARGET_GENE_SYMBOL'] #4779 SOMAmers targeting 4102 unique proteins
colnames(proteomics_human) <- name_symbol

### PCA for all probes of proteins
pc_protein <- read.csv("Proteomics_PC.tsv", sep = '\t') %>%
  mutate(PATNO = paste0("PP-", PATNO))
# pc_protein <- prcomp(t(proteomics_human), scale = T)
# pc4_protein <- pc_protein$rotation[, 1:4]
# pc4_protein <- cbind("PATNO"=proteomics_df[, 1], pc4_protein)

# merge PCs with phenotype
df.pheno <- merge(pc_protein, df_pheno, by.x = "PATNO", by.y = "row.names") %>%
  column_to_rownames("PATNO") #373

# extract common IDs in proteomics and phenotype
df_proteomics <- proteomics_human[rownames(proteomics_human) %in% rownames(df.pheno),]

mod <- model.matrix(~ cluster-1 + AGE_AT_VISIT + SEX + duration_yrs + type_mutation + PC1 + PC2 + PC3 + PC4, data = df.pheno)
colnames(mod)[1:7] <- c('C1', 'C2', 'AGE_AT_VISIT', 'Male', 'duration_yrs', 'GBA1', 'LRRK2')
res <- de_func(model = mod, feature = df_proteomics, contrasts = "C1-C2", 
               table.out = "/Volumes/Extreme SSD/work_bidmc/PPMI/iCluster/PD/DE/DE_icluster_PD_proteomics_csf.tsv")
# no significant CSF proteins

###############################################################################
#                           End: CSF Proteomics                               #
###############################################################################
