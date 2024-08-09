library(tidyverse)
library(icesTAF)

rm(list = ls())
gc()

main_dir <- "C:/Users/vasileioubill95/Desktop/Pipeline"
mkdir(sprintf("%s/Masterfiles/",main_dir))

## Activity Master file - Case Study 1 ##

# MGH #

mgh <- read.csv(sprintf("%s/Case_study_1/5.Results/Activity/MGH/DA.csv", main_dir), header = T)
mgh_DE <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/MGH/DE.csv", main_dir), header = T)
mgh <- mgh[order(mgh$ID ,decreasing=FALSE),]

mgh <- mgh %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "TRUE", "Positive")
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "FALSE", "Negative")

mgh_DE_sig <- filter(mgh_DE,adj.P.Val < 0.05)
mgh_overt <- mgh[mgh$ID %in% mgh_DE_sig$ID,]
mgh_overt$Type <- "Overt"
mgh_hidden <- mgh[!(mgh$ID %in% mgh_DE_sig$ID),]
mgh_hidden$Type <- "Hidden"

mgh <- rbind(mgh_overt, mgh_hidden)

mgh <- mgh %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mgh$Significant <- str_replace_all(mgh$Significant, "TRUE", "Sig")
mgh$Significant <- str_replace_all(mgh$Significant, "FALSE", "NonSig")

mgh$Type[mgh$Significant == "NonSig"] <- "NA"
rownames(mgh) <- mgh$ID
mgh <- subset(mgh, select = c(adj.P.Val, Z.statistics, Perturbational_direct, Type))
colnames(mgh) <- c("adj.P.Val_MGH", "Z.statistics_MGH", "Perturbational_direct_MGH", "Type_MGH")

# Mayo #

mayo <- read.csv(sprintf("%s/Case_study_1/5.Results/Activity/Mayo/DA.csv", main_dir), header = T)
mayo_DE <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/Mayo/DE.csv", main_dir), header = T)
mayo <- mayo[order(mayo$ID ,decreasing=FALSE),]

mayo <- mayo %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mayo$Perturbational_direct <- str_replace_all(mayo$Perturbational_direct, "TRUE", "Positive")
mayo$Perturbational_direct <- str_replace_all(mayo$Perturbational_direct, "FALSE", "Negative")

mayo_DE_sig <- filter(mayo_DE,adj.P.Val < 0.05)
mayo_overt <- mayo[mayo$ID %in% mayo_DE_sig$ID,]
mayo_overt$Type <- "Overt"
mayo_hidden <- mayo[!(mayo$ID %in% mayo_DE_sig$ID),]
mayo_hidden$Type <- "Hidden"

mayo <- rbind(mayo_overt, mayo_hidden)

mayo <- mayo %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mayo$Significant <- str_replace_all(mayo$Significant, "TRUE", "Sig")
mayo$Significant <- str_replace_all(mayo$Significant, "FALSE", "NonSig")

mayo$Type[mayo$Significant == "NonSig"] <- "NA"
rownames(mayo) <- mayo$ID
mayo <- subset(mayo, select = c(adj.P.Val, Z.statistics, Perturbational_direct, Type))
colnames(mayo) <- c("adj.P.Val_Mayo", "Z.statistics_Mayo", "Perturbational_direct_Mayo", "Type_Mayo")

# Stanford #

stanford <- read.csv(sprintf("%s/Case_study_1/5.Results/Activity/Stanford/DA.csv", main_dir), header = T)
stanford_DE <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/Stanford/DE.csv", main_dir), header = T)
stanford <- stanford[order(stanford$ID ,decreasing=FALSE),]

stanford <- stanford %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
stanford$Perturbational_direct <- str_replace_all(stanford$Perturbational_direct, "TRUE", "Positive")
stanford$Perturbational_direct <- str_replace_all(stanford$Perturbational_direct, "FALSE", "Negative")

stanford_DE_sig <- filter(stanford_DE,adj.P.Val < 0.05)
stanford_overt <- stanford[stanford$ID %in% stanford_DE_sig$ID,]
stanford_overt$Type <- "Overt"
stanford_hidden <- stanford[!(stanford$ID %in% stanford_DE_sig$ID),]
stanford_hidden$Type <- "Hidden"

stanford <- rbind(stanford_overt, stanford_hidden)

stanford <- stanford %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
stanford$Significant <- str_replace_all(stanford$Significant, "TRUE", "Sig")
stanford$Significant <- str_replace_all(stanford$Significant, "FALSE", "NonSig")

stanford$Type[stanford$Significant == "NonSig"] <- "NA"
rownames(stanford) <- stanford$ID
stanford <- subset(stanford, select = c(adj.P.Val, Z.statistics, Perturbational_direct, Type))
colnames(stanford) <- c("adj.P.Val_Stanford", "Z.statistics_Stanford", "Perturbational_direct_Stanford", "Type_Stanford")

rm(list=setdiff(ls(), c("main_dir", "mgh", "mayo", "stanford")))

# Merge #

mgh_mayo <- merge(mgh, mayo, by = 0)
rownames(mgh_mayo) <- mgh_mayo$Row.names
mgh_mayo <- subset(mgh_mayo, select = -c(Row.names))

all_df <- merge(mgh_mayo, stanford, by = 0)
rownames(all_df) <- all_df$Row.names
all_df <- subset(all_df , select = -c(Row.names))


cluster <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_OLINK_NPX.txt", main_dir), header = T, sep = ";")
cluster <- subset(cluster, select = c(Assay,Panel))
cluster <- cluster[!duplicated(cluster$Assay), ]
rownames(cluster) <- cluster$Assay
cluster <- subset(cluster, select = -c(Assay))

df_cluster <- merge(all_df, cluster, by = 0)
rownames(df_cluster) <- df_cluster$Row.names
names(df_cluster)[names(df_cluster) == 'Row.names'] <- 'ID'

APNet_masterfile <- df_cluster

# merge Hidden

APNet_masterfile_yes <- filter(APNet_masterfile, Type_MGH == "Hidden" | Type_Mayo == "Hidden" | Type_Stanford == "Hidden")
APNet_masterfile_no <- APNet_masterfile[!(APNet_masterfile$ID %in% APNet_masterfile_yes$ID),]
APNet_masterfile_nosig <- filter(APNet_masterfile_no, Type_MGH == "NA" & Type_Mayo == "NA" & Type_Stanford == "NA")
APNet_masterfile_no_sig <- APNet_masterfile_no[!(APNet_masterfile_no$ID %in% APNet_masterfile_nosig$ID),]

APNet_masterfile_yes$Total_Hidden <- "Hidden"
APNet_masterfile_no_sig$Total_Hidden <- "Overt"
APNet_masterfile_nosig$Total_Hidden <- "NA"

APNet_masterfile <- rbind(APNet_masterfile_yes, APNet_masterfile_no_sig, APNet_masterfile_nosig)

rm(list=setdiff(ls(), c("main_dir","APNet_masterfile")))

APNet_masterfile$Perturbational_direct_Total <- "NA"

APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Positive" & APNet_masterfile$Perturbational_direct_Mayo == "Positive" & APNet_masterfile$Perturbational_direct_Stanford == "Positive"] <- "Positive"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Negative" & APNet_masterfile$Perturbational_direct_Mayo == "Negative" & APNet_masterfile$Perturbational_direct_Stanford == "Negative"] <- "Negative"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_Total == "NA"] <- "Mixed"

common <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Overlaps_Activity/common.txt", main_dir), header = T)
common <- common$x

APNet_masterfile <- APNet_masterfile[APNet_masterfile$ID %in% common,]
# save

write.csv(APNet_masterfile, file = sprintf("%s/Masterfiles/APNet_masterfile.csv", main_dir), row.names = F, quote = F)

## Expression Master file - Case Study 1 ##

rm(list=setdiff(ls(), "main_dir"))
gc()

# MGH #

mgh <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/MGH/DE.csv", main_dir), header = T)
mgh <- mgh[order(mgh$ID ,decreasing=FALSE),]

mgh <- mgh %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "TRUE", "Positive")
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "FALSE", "Negative")

mgh <- mgh %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mgh$Significant <- str_replace_all(mgh$Significant, "TRUE", "Sig")
mgh$Significant <- str_replace_all(mgh$Significant, "FALSE", "NonSig")

mgh$Type[mgh$Significant == "NonSig"] <- "NA"
rownames(mgh) <- mgh$ID
mgh <- subset(mgh, select = c(adj.P.Val, Z.statistics, Perturbational_direct))
colnames(mgh) <- c("adj.P.Val_MGH", "Z.statistics_MGH", "Perturbational_direct_MGH")

# Mayo #

mayo <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/Mayo/DE.csv", main_dir), header = T)
mayo <- mayo[order(mayo$ID ,decreasing=FALSE),]

mayo <- mayo %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mayo$Perturbational_direct <- str_replace_all(mayo$Perturbational_direct, "TRUE", "Positive")
mayo$Perturbational_direct <- str_replace_all(mayo$Perturbational_direct, "FALSE", "Negative")

mayo <- mayo %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mayo$Significant <- str_replace_all(mayo$Significant, "TRUE", "Sig")
mayo$Significant <- str_replace_all(mayo$Significant, "FALSE", "NonSig")

mayo$Type[mayo$Significant == "NonSig"] <- "NA"
rownames(mayo) <- mayo$ID
mayo <- subset(mayo, select = c(adj.P.Val, Z.statistics, Perturbational_direct))
colnames(mayo) <- c("adj.P.Val_Mayo", "Z.statistics_Mayo", "Perturbational_direct_Mayo")

# Stanford #

stanford <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/Stanford/DE.csv", main_dir), header = T)
stanford <- stanford[order(stanford$ID ,decreasing=FALSE),]

stanford <- stanford %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
stanford$Perturbational_direct <- str_replace_all(stanford$Perturbational_direct, "TRUE", "Positive")
stanford$Perturbational_direct <- str_replace_all(stanford$Perturbational_direct, "FALSE", "Negative")

stanford <- stanford %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
stanford$Significant <- str_replace_all(stanford$Significant, "TRUE", "Sig")
stanford$Significant <- str_replace_all(stanford$Significant, "FALSE", "NonSig")

stanford$Type[stanford$Significant == "NonSig"] <- "NA"
rownames(stanford) <- stanford$ID
stanford <- subset(stanford, select = c(adj.P.Val, Z.statistics, Perturbational_direct))
colnames(stanford) <- c("adj.P.Val_Standord", "Z.statistics_Stanford", "Perturbational_direct_Stanford")

# Merged #

mgh_mayo <- merge(mgh, mayo, by = 0)
rownames(mgh_mayo) <- mgh_mayo$Row.names
mgh_mayo <- subset(mgh_mayo, select = -c(Row.names))

all_df <- merge(mgh_mayo, stanford, by = 0)
rownames(all_df) <- all_df$Row.names
all_df <- subset(all_df , select = -c(Row.names))

cluster <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_OLINK_NPX.txt", main_dir), header = T, sep = ";")
cluster <- subset(cluster, select = c(Assay,Panel))
cluster <- cluster[!duplicated(cluster$Assay), ]
rownames(cluster) <- cluster$Assay
cluster <- subset(cluster, select = -c(Assay))

df_cluster <- merge(all_df, cluster, by = 0)
rownames(df_cluster) <- df_cluster$Row.names
names(df_cluster)[names(df_cluster) == 'Row.names'] <- 'ID'

APNet_masterfile <- df_cluster

rm(list=setdiff(ls(), c("main_dir","APNet_masterfile")))

APNet_masterfile$Perturbational_direct_Total <- "NA"

APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Positive" & APNet_masterfile$Perturbational_direct_Mayo == "Positive" & APNet_masterfile$Perturbational_direct_Stanford == "Positive"] <- "Positive"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Negative" & APNet_masterfile$Perturbational_direct_Mayo == "Negative" & APNet_masterfile$Perturbational_direct_Stanford == "Negative"] <- "Negative"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_Total == "NA"] <- "Mixed"

common <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Overlaps_Expression/common.txt", main_dir), header = T)
common <- common$x

APNet_masterfile <- APNet_masterfile[APNet_masterfile$ID %in% common,]

# save

write.csv(APNet_masterfile, file = sprintf("%s/Masterfiles/APNet_masterfile_DE.csv", main_dir), row.names = F, quote = F)

## Activity Master file - Case Study 2 ##

# MGH #

mgh <- read.csv(sprintf("%s/Case_study_1/5.Results/Activity/MGH/DA.csv", main_dir), header = T)
mgh_DE <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/MGH/DE.csv", main_dir), header = T)
mgh <- mgh[order(mgh$ID ,decreasing=FALSE),]

mgh <- mgh %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "TRUE", "Positive")
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "FALSE", "Negative")

mgh_DE_sig <- filter(mgh_DE,adj.P.Val < 0.05)
mgh_overt <- mgh[mgh$ID %in% mgh_DE_sig$ID,]
mgh_overt$Type <- "Overt"
mgh_hidden <- mgh[!(mgh$ID %in% mgh_DE_sig$ID),]
mgh_hidden$Type <- "Hidden"

mgh <- rbind(mgh_overt, mgh_hidden)

mgh <- mgh %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mgh$Significant <- str_replace_all(mgh$Significant, "TRUE", "Sig")
mgh$Significant <- str_replace_all(mgh$Significant, "FALSE", "NonSig")

mgh$Type[mgh$Significant == "NonSig"] <- "NA"
rownames(mgh) <- mgh$ID
mgh <- subset(mgh, select = c(adj.P.Val, Z.statistics, Perturbational_direct, Type))
colnames(mgh) <- c("adj.P.Val_MGH", "Z.statistics_MGH", "Perturbational_direct_MGH", "Type_MGH")

# scMGH #

scmgh <- read.csv(sprintf("%s/Case_study_2/4.Results/DAG_tf.csv", main_dir), header = T)
scmgh_DE <- read.csv(sprintf("%s/Case_study_2/4.Results/DEP.csv", main_dir), header = T)
scmgh <- scmgh[order(scmgh$id ,decreasing=FALSE),]

scmgh <- scmgh %>% mutate(Perturbational_direct=ifelse(Z_Severe > 0,T,F))
scmgh$Perturbational_direct <- str_replace_all(scmgh$Perturbational_direct, "TRUE", "Positive")
scmgh$Perturbational_direct <- str_replace_all(scmgh$Perturbational_direct, "FALSE", "Negative")

scmgh$id <- str_remove(scmgh$id, ".TF")

common <- read.table(sprintf("%s/Case_study_2/5.Overlaps/Overlaps/common.txt", main_dir), header = T)
common <- common$x

scmgh <- scmgh[scmgh$id %in% common,]
scmgh_DE <- scmgh_DE[scmgh_DE$id %in% common,]

scmgh_DE_sig <- filter(scmgh_DE,pval_Severe < 0.05)
scmgh_overt <- scmgh[scmgh$id %in% scmgh_DE_sig$id,]
scmgh_overt$Type <- "Overt"
scmgh_hidden <- scmgh[!(scmgh$id %in% scmgh_DE_sig$id),]
scmgh_hidden$Type <- "Hidden"

scmgh <- rbind(scmgh_overt, scmgh_hidden)

scmgh <- scmgh %>% mutate(Significant=ifelse(pval_Severe < 0.05,T,F))
scmgh$Significant <- str_replace_all(scmgh$Significant, "TRUE", "Sig")
scmgh$Significant <- str_replace_all(scmgh$Significant, "FALSE", "NonSig")

scmgh$Type[scmgh$Significant == "NonSig"] <- "NA"
rownames(scmgh) <- scmgh$id
scmgh <- subset(scmgh, select = c(pval_Severe, Z_Severe, Perturbational_direct, Type))
colnames(scmgh) <- c("adj.P.Val_scMGH", "Z.statistics_scMGH", "Perturbational_direct_scMGH", "Type_scMGH")

rm(list=setdiff(ls(), c("scmgh", "mgh", "main_dir")))

# Merged #

all_df <- merge(mgh, scmgh, by = 0)
rownames(all_df) <- all_df$Row.names
all_df <- subset(all_df , select = -c(Row.names))

cluster <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_OLINK_NPX.txt", main_dir), header = T, sep = ";")
cluster <- subset(cluster, select = c(Assay,Panel))
cluster <- cluster[!duplicated(cluster$Assay), ]
rownames(cluster) <- cluster$Assay
cluster <- subset(cluster, select = -c(Assay))

df_cluster <- merge(all_df, cluster, by = 0)
rownames(df_cluster) <- df_cluster$Row.names
names(df_cluster)[names(df_cluster) == 'Row.names'] <- 'ID'

APNet_masterfile <- df_cluster

rm(list=setdiff(ls(), c("main_dir","APNet_masterfile")))

APNet_masterfile$Perturbational_direct_Total <- "NA"

APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Positive" & APNet_masterfile$Perturbational_direct_Mayo == "Positive" & APNet_masterfile$Perturbational_direct_Stanford == "Positive"] <- "Positive"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Negative" & APNet_masterfile$Perturbational_direct_Mayo == "Negative" & APNet_masterfile$Perturbational_direct_Stanford == "Negative"] <- "Negative"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_Total == "NA"] <- "Mixed"

write.csv(APNet_masterfile, file = sprintf("%s/Masterfiles/APNet_masterfile_scMGH.csv", main_dir), row.names = F, quote = F)

## Expression Master file - Case Study 2 ##

rm(list=setdiff(ls(), "main_dir"))
gc()

# MGH #

mgh <- read.csv(sprintf("%s/Case_study_1/5.Results/Expression/MGH/DE.csv", main_dir), header = T)
mgh <- mgh[order(mgh$ID ,decreasing=FALSE),]

mgh <- mgh %>% mutate(Perturbational_direct=ifelse(Z.statistics > 0,T,F))
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "TRUE", "Positive")
mgh$Perturbational_direct <- str_replace_all(mgh$Perturbational_direct, "FALSE", "Negative")

mgh <- mgh %>% mutate(Significant=ifelse(adj.P.Val < 0.05,T,F))
mgh$Significant <- str_replace_all(mgh$Significant, "TRUE", "Sig")
mgh$Significant <- str_replace_all(mgh$Significant, "FALSE", "NonSig")

mgh$Type[mgh$Significant == "NonSig"] <- "NA"
rownames(mgh) <- mgh$ID
mgh <- subset(mgh, select = c(adj.P.Val, Z.statistics, Perturbational_direct))
colnames(mgh) <- c("adj.P.Val_MGH", "Z.statistics_MGH", "Perturbational_direct_MGH")

# scMGH #

scmgh <- read.csv(sprintf("%s/Case_study_2/4.Results/DEP.csv", main_dir), header = T)
scmgh <- scmgh[order(scmgh$id ,decreasing=FALSE),]

scmgh <- scmgh %>% mutate(Perturbational_direct=ifelse(Z_Severe > 0,T,F))
scmgh$Perturbational_direct <- str_replace_all(scmgh$Perturbational_direct, "TRUE", "Positive")
scmgh$Perturbational_direct <- str_replace_all(scmgh$Perturbational_direct, "FALSE", "Negative")

scmgh$id <- str_remove(scmgh$id, ".TF")

common <- read.table(sprintf("%s/Case_study_2/5.Overlaps/Overlaps/common_DE.txt", main_dir), header = T)
common <- common$x

scmgh <- scmgh[scmgh$id %in% common,]

scmgh <- scmgh %>% mutate(Significant=ifelse(pval_Severe < 0.05,T,F))
scmgh$Significant <- str_replace_all(scmgh$Significant, "TRUE", "Sig")
scmgh$Significant <- str_replace_all(scmgh$Significant, "FALSE", "NonSig")

scmgh$Type[scmgh$Significant == "NonSig"] <- "NA"
rownames(scmgh) <- scmgh$id
scmgh <- subset(scmgh, select = c(pval_Severe, Z_Severe, Perturbational_direct, Type))
colnames(scmgh) <- c("adj.P.Val_scMGH", "Z.statistics_scMGH", "Perturbational_direct_scMGH", "Type_scMGH")

rm(list=setdiff(ls(), c("scmgh", "mgh", "main_dir")))
# Merged #

all_df <- merge(mgh, scmgh, by = 0)
rownames(all_df) <- all_df$Row.names
all_df <- subset(all_df , select = -c(Row.names))

cluster <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_OLINK_NPX.txt", main_dir), header = T, sep = ";")
cluster <- subset(cluster, select = c(Assay,Panel))
cluster <- cluster[!duplicated(cluster$Assay), ]
rownames(cluster) <- cluster$Assay
cluster <- subset(cluster, select = -c(Assay))

df_cluster <- merge(all_df, cluster, by = 0)
rownames(df_cluster) <- df_cluster$Row.names
names(df_cluster)[names(df_cluster) == 'Row.names'] <- 'ID'

APNet_masterfile <- df_cluster

rm(list=setdiff(ls(), c("main_dir","APNet_masterfile")))

APNet_masterfile$Perturbational_direct_Total <- "NA"

APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Positive" & APNet_masterfile$Perturbational_direct_Mayo == "Positive" & APNet_masterfile$Perturbational_direct_Stanford == "Positive"] <- "Positive"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_MGH == "Negative" & APNet_masterfile$Perturbational_direct_Mayo == "Negative" & APNet_masterfile$Perturbational_direct_Stanford == "Negative"] <- "Negative"
APNet_masterfile$Perturbational_direct_Total[APNet_masterfile$Perturbational_direct_Total == "NA"] <- "Mixed"

# save

write.csv(APNet_masterfile, file = sprintf("%s/Masterfiles/APNet_masterfile_scMGH_DE.csv", main_dir), row.names = F, quote = F)


