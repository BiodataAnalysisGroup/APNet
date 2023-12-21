#Script for Figures in APNet

# load libraries

library(tidyverse)
library(ggvenn)
library(patchwork)
library(icesTAF)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(openxlsx)


# create folders per Figure
main_dir <- "C:/Users/vasileioubill95/Desktop/Pipeline"
mkdir("Figures")

## FIGURE 3 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Figure_3", main_dir))

# create VennDiagrams from common files

mgh_neg <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/MGH/MGH_neg.txt", main_dir), header = T)
mgh_neg <- mgh_neg$x
mgh_pos <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/MGH/MGH_pos.txt", main_dir), header = T)
mgh_pos <- mgh_pos$x

mayo_neg <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/Mayo/Mayo_neg.txt", main_dir), header = T)
mayo_neg <- mayo_neg$x
mayo_pos <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/Mayo/Mayo_pos.txt", main_dir), header = T)
mayo_pos <- mayo_pos$x

stanford_neg <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/Stanford/Stanford_neg.txt", main_dir), header = T)
stanford_neg <- stanford_neg$x
stanford_pos <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/Stanford/Stanford_pos.txt", main_dir), header = T)
stanford_pos <- stanford_pos$x

pos <- list(MGH=sample(mgh_pos),
             Mayo=sample(mayo_pos),
             Stanford=sample(stanford_pos))

neg <- list(MGH=sample(mgh_neg),
             Mayo=sample(mayo_neg),
             Stanford=sample(stanford_neg))

p_pos <- ggvenn::ggvenn(pos, set_name_size = 4, text_size = 4, show_percentage = FALSE) + 
  labs(title = "Common Positive DA") +
  theme(text = element_text(size=10))

p_neg <-  ggvenn::ggvenn(neg, set_name_size = 4, text_size = 4, show_percentage = FALSE) + 
  labs(title = "Common Negative DA") +
  theme(text = element_text(size=10))

p_venn <- p_pos + p_neg

ggsave(sprintf("%s/Figures/Figure_3/venndiagram.svg",main_dir), plot = p_venn, device = "svg", width = 8, height = 7)

# create Bubble plot from EnrichR KG downloaded for Positive and Negative Drivers

# Positive #

enrich_pos <- read.csv(sprintf("%s/Case_study_1/EnrichR_KG/Activity/Enrichr-KG_pos.csv",main_dir))
enrich_pos <- enrich_pos[order(enrich_pos$q.value, decreasing = FALSE), ]

enrich_pos$Term <- str_remove(enrich_pos$Term, ".WP.*")
enrich_pos$Term <- str_remove(enrich_pos$Term, "..GO.*")
enrich_pos$Term <- str_remove(enrich_pos$Term, ".R.HSA.*")
enrich_pos$Term <- str_replace_all(enrich_pos$Term, "\\.", " ")

enrich_pos <- enrich_pos[!duplicated(enrich_pos[,c("Term")]),]

enrich_pos <- enrich_pos[1:30,]

enrich_pos$logqvalue <- -log10(enrich_pos$q.value)
enrich_pos <- enrich_pos[order(enrich_pos$logqvalue, decreasing = F),]

enrich_pos$Term[[11]] <- "SARS-CoV-2 innate immunity evasion"

# Plot
p_bubble_pos <- ggplot(data = enrich_pos, aes(x=logqvalue, y=Term, size = combined.score , color = z.score )) +
  geom_point(alpha=1) +
  ggtitle("Pathways DA Up-regulated") +
  ylab("Terms") +
  xlab("-10log(P Value)") +
  labs(color = "Z Score") +
  labs(size = "Combine Score") +
  scale_y_discrete(limits = enrich_pos$Term) +
  scale_colour_continuous(trans = 'reverse') +
  scale_color_gradient(low = "firebrick1", high = "firebrick")

ggsave(sprintf("%s/Figures/Figure_3/bubble_plot_pos.svg",main_dir), plot = p_bubble_pos, device = "svg", width = 8, height = 4.5)

# Negative #

enrich_neg <- read.csv(sprintf("%s/Case_study_1/EnrichR_KG/Activity/Enrichr-KG_neg.csv",main_dir))
enrich_neg <- enrich_neg[order(enrich_neg$q.value, decreasing = FALSE), ]

enrich_neg$Term <- str_remove(enrich_neg$Term, ".WP.*")
enrich_neg$Term <- str_remove(enrich_neg$Term, "..GO.*")
enrich_neg$Term <- str_remove(enrich_neg$Term, ".R.HSA.*")
enrich_neg$Term <- str_replace_all(enrich_neg$Term, "\\.", " ")

enrich_neg <- enrich_neg[!duplicated(enrich_neg[,c("Term")]),]

enrich_neg <- enrich_neg[1:30,]

enrich_neg$logqvalue <- -log10(enrich_neg$q.value)
enrich_neg <- enrich_neg[order(enrich_neg$logqvalue, decreasing = F),]

enrich_neg$Term[[16]] <- "homophilic cell adhesion"

# Plot
p_bubble_neg <- ggplot(data = enrich_neg, aes(x=logqvalue, y=Term, size = combined.score , color = z.score )) +
  geom_point(alpha=1) +
  ggtitle("Pathways DA Up-regulated") +
  ylab("Terms") +
  xlab("-10log(P Value)") +
  labs(color = "Z Score") +
  labs(size = "Combine Score") +
  scale_y_discrete(limits = enrich_neg$Term) +
  scale_colour_continuous(trans = 'reverse') +
  scale_color_gradient(low = "royalblue1", high = "royalblue4")

ggsave(sprintf("%s/Figures/Figure_3/bubble_plot_neg.svg", main_dir), plot = p_bubble_neg, device = "svg", width = 8, height = 4.5)

## FIGURE 4 ##

## RUN first the auc_plot.py

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Figure_4", main_dir))

# create auc from APNet approach

# Mayo #

df <- read.table(sprintf("%s/Figures/tables_auc/df_auc.csv", main_dir),
                   header = T, sep = ",")

mayo_auc <- filter(df, Model == "APNET_Mayo (AUC: 0.96)")
 
p_mayo_auc <- ggplot(mayo_auc, aes(x = FPR, y = TPR, fill = Model, color = Model)) +
  geom_line(aes(color=Model)) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate', title = 'MAYO AUC: 0.96 - F1: 0.91') +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(sprintf("%s/Figures/Figure_4/mayo_auc.svg", main_dir), plot = p_mayo_auc, device = "svg", width = 5, height = 5)

# Staford #

df <- read.table(sprintf("%s/Figures/tables_auc/df_auc.csv",main_dir),
                 header = T, sep = ",")

stanford_auc <- filter(df, Model == "APNet_Stanford (AUC: 0.91)")


p_stanford_auc <- ggplot(stanford_auc, aes(x = FPR, y = TPR, fill = Model, color = Model)) +
  geom_line(aes(color=Model)) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate', title = 'STANFORD AUC: 0.91 - F1: 0.68') +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(sprintf("%s/Figures/Figure_4/stanford_auc.svg",main_dir), plot = p_stanford_auc, device = "svg", width = 5, height = 5)

# create heatmaps from top20 SHAP values and pathways sc1 APNet output

# Mayo #

mayo_ac <- as.vector(c("JUN", "IL6", "MAPK9", "LYN", "TNFRSF1A", "AREG",
                       "NTF4", "NCF2", "TNFRSF10A", "HGF", "FLT3", "CKAP4",
                       "FLT3LG", "SDC1", "EFNA1", "TNFRSF10B", "TNFSF11",
                       "ACAA1", "TIMP1", "PTEN"))


df <- read.xlsx(sprintf("%s/Case_study_1/7.PASNet/Activity/Output/Mayo_test/sc1_weights_fixed.xlsx",main_dir), rowNames = TRUE)
df <- df[,colnames(df) %in% mayo_ac]

df_1 <- as.data.frame(abs(df))
df_1 <- df_1[order(rowSums(df_1),decreasing=T),]
df_1 <- df_1[,order(colSums(df_1), decreasing = T)]
df_1$names <- rownames(df_1)
df_1 <- filter(df_1, names != "Apoptosis WP254")
df_1 <- subset(df_1, select = -c(names))

df_1 <- df_1[1:50,]
vector <- rownames(df_1)

df <- df[rownames(df) %in% vector,]

df <- df[order(rowSums(df),decreasing=T),]
df <- df[,order(colSums(df), decreasing = T)]

rownames(df) <- str_remove(rownames(df), ".WP.*")
rownames(df) <- str_remove(rownames(df), "..GO.*")
rownames(df) <- str_remove(rownames(df), ".R.HSA.*")
rownames(df) <- str_replace_all(rownames(df), "\\.", " ")

rownames(df)[[24]] <- "Caspase Activation Via Death Receptor"
rownames(df)[[30]] <- "HSP90-Chaperone Cycle for SHR"

df <- as.matrix(df)
col_fun = colorRamp2(c(min(df), 0, max(df)), c("blue", "white", "red"))
col_fun(seq(min(df), max(df)))

p_mayo_heatmap <- as.ggplot(Heatmap(df, 
                               name = "Weights",
                               column_title = "Proteins", row_title = "Pathways", cluster_columns = F,
                               cluster_rows = F, col = col_fun,
                               row_names_gp = gpar(fontsize = 7.5), 
                               column_names_gp = gpar(fontsize = 9)
))

ggsave(sprintf("%s/Figures/Figure_4/sc1_heatmap_mayo_ac_top20_SHAP.svg",main_dir), plot = p_mayo_heatmap, device = "svg", width = 8, height = 7)

# Stanford #

stanford_ac <- as.vector(c("PTEN", "TNFRSF1A", "IL6", "BAX", "LYN", "LTA", 
                           "KDR", "COL1A1", "JUN", "CCL7", "TNFRSF10B", 
                           "TNFRSF10A", "TNFSF10", "EGFR", "ERBB2", "CCL22",
                           "PODXL", "SEMA4D", "KIT", "ROBO1"))

df <- read.xlsx(sprintf("%s/Case_study_1/7.PASNet/Activity/Output/Stanford_test/sc1_weights.xlsx",main_dir), rowNames = TRUE)
df <- df[,colnames(df) %in% stanford_ac]

df_1 <- as.data.frame(abs(df))
df_1 <- df_1[order(rowSums(df_1),decreasing=T),]
df_1 <- df_1[,order(colSums(df_1), decreasing = T)]

df_1 <- df_1[1:50,]
vector <- rownames(df_1)

df <- df[rownames(df) %in% vector,]

df <- df[order(rowSums(df),decreasing=T),]
df <- df[,order(colSums(df), decreasing = T)]

rownames(df) <- str_remove(rownames(df), ".WP.*")
rownames(df) <- str_remove(rownames(df), "..GO.*")
rownames(df) <- str_remove(rownames(df), ".R.HSA.*")
rownames(df) <- str_replace_all(rownames(df), "\\.", " ")

rownames(df)[35] <- "TWEAK Signaling Pathway"
rownames(df)[25] <- "Photodynamic therapy-induced AP-1 survival signal"
  
df <- as.matrix(df)
col_fun = colorRamp2(c(min(df), 0, max(df)), c("blue", "white", "red"))
col_fun(seq(min(df), max(df)))

p_stanford_heatmap <- as.ggplot(Heatmap(df, 
                               name = "Weights",
                               column_title = "Proteins", row_title = "Pathways", cluster_columns = F,
                               cluster_rows = F, col = col_fun,
                               row_names_gp = gpar(fontsize = 7.5), 
                               column_names_gp = gpar(fontsize = 9)
))

ggsave(sprintf("%s/Figures/Figure_4/sc1_heatmap_stanford_ac_top20_SHAP.svg",main_dir), plot = p_stanford_heatmap, device = "svg", width = 8, height = 7)

# create heatmap with top20 SHAP values from PASNet by using MGH activity matrix

# Mayo SHAP #

vector <- as.vector(c("JUN", "IL6", "MAPK9", "LYN", "TNFRSF1A", "AREG",
                       "NTF4", "NCF2", "TNFRSF10A", "HGF", "FLT3", "CKAP4",
                       "FLT3LG", "SDC1", "EFNA1", "TNFRSF10B", "TNFSF11",
                       "ACAA1", "TIMP1", "PTEN"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH APNet - Mayo top20 SHAP values")


ggsave(sprintf("%s/Figures/Figure_4/mgh_mayo_heatmap_activity.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford SHAP #

vector <- as.vector(c("PTEN", "TNFRSF1A", "IL6", "BAX", "LYN", "LTA", 
                           "KDR", "COL1A1", "JUN", "CCL7", "TNFRSF10B", 
                           "TNFRSF10A", "TNFSF10", "EGFR", "ERBB2", "CCL22",
                           "PODXL", "SEMA4D", "KIT", "ROBO1"))


df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH APNet - Stanford top20 SHAP values")

ggsave(sprintf("%s/Figures/Figure_4/mgh_stanford_heatmap_activity.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# create heatmap with top20 SHAP values by using Mayo and Stanford activity matrices

# Mayo #

vector <- as.vector(c("JUN", "IL6", "MAPK9", "LYN", "TNFRSF1A", "AREG",
                      "NTF4", "NCF2", "TNFRSF10A", "HGF", "FLT3", "CKAP4",
                      "FLT3LG", "SDC1", "EFNA1", "TNFRSF10B", "TNFSF11",
                      "ACAA1", "TIMP1", "PTEN"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Mayo/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]

pheno <- read.xlsx(sprintf("%s/Case_study_1/1.Raw_files/Mayo/MC_raw.xlsx",main_dir))
rownames(pheno) <- pheno$Sample

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]

df <- df[,colnames(df) %in% rownames(pheno)]

matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- subset(pheno, select = c(WHOscale))
colnames(pheno)[colnames(pheno) == "WHOscale"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "1", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "2", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "3", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "4", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "5", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "6", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "7", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "8", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(Severe="red", NonSevere="blue"))
)

p <- Heatmap(matrix_norm, km = 1, 
                  show_row_names = T, 
                  show_column_names = F, 
                  cluster_rows = T, 
                  cluster_columns = T, 
                  name = "Activity", 
                  top_annotation = ha,
                  column_title = "Mayo APNet - Mayo top20 SHAP values")


ggsave(sprintf("%s/Figures/Figure_4/mayo_heatmap_activity.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford #

vector <- as.vector(c("PTEN", "TNFRSF1A", "IL6", "BAX", "LYN", "LTA", 
                      "KDR", "COL1A1", "JUN", "CCL7", "TNFRSF10B", 
                      "TNFRSF10A", "TNFSF10", "EGFR", "ERBB2", "CCL22",
                      "PODXL", "SEMA4D", "KIT", "ROBO1"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Stanford/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")

pheno <- read.csv(sprintf("%s/Case_study_1/1.Raw_files/Stanford/PatientCharacteristics.csv",main_dir),
                  header = T)
rownames(pheno) <- pheno$sampleID

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]

matrix <- as.matrix(df)

matrix_norm <- scale(t(matrix), center = T, scale = T)
matrix_norm <- matrix_norm[order(rownames(matrix_norm)),]
matrix_norm <- t(matrix_norm)

pheno <- pheno[order(rownames(pheno)),]
pheno <- subset(pheno, select = c(Severity))
colnames(pheno)[colnames(pheno) == "Severity"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "Asymptomatic", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Healthy", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Mild", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Moderate", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Severe", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(NonSevere="blue", Severe="red"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "Stanford APNet - Stanford top20 SHAP values"
             )

ggsave(sprintf("%s/Figures/Figure_4/stanford_heatmap_activity.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

## FIGURE 5 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Figure_5", main_dir))

# create auc plot from APNet approach

# scMGH #

df <- read.table(sprintf("%s/Figures/tables_auc/df_auc.csv",main_dir),
                 header = T, sep = ",")

scmgh_auc <- filter(df, Model == "APNet_scMGH (AUC: 0.99)")

p_scmgh_auc <- ggplot(scmgh_auc, aes(x = FPR, y = TPR, fill = Model, color = Model)) +
  geom_line(aes(color=Model)) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate', title = 'scMGH AUC: 0.99 - F1: 0.975') +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(sprintf("%s/Figures/Figure_5/scmgh_auc.svg", main_dir), plot = p_scmgh_auc, device = "svg", width = 5, height = 5)

# create heatmaps from top20 SHAP values and pathways sc1 APNet output

# scMGH #

vector <- c("JUN", "TIMP1", "S100A11", "ITGB2","ITGAM", "IKBKG", "S100A12", "IL6", 
               "CD63", "LAMP2", "BIRC2", "HMOX1", "LGALS1", "NFATC1", "IL10RA", "ATP6AP2",
               "CD4", "TNFSF10", "MAPK9", "ITGB1")


df <- read.xlsx(sprintf("%s/Case_study_2/7.PASNet/Output/sc1_weights.xlsx", main_dir), rowNames = TRUE)
df <- df[,colnames(df) %in% vector]

df_1 <- as.data.frame(abs(df))
df_1 <- df_1[order(rowSums(df_1),decreasing=T),]
df_1 <- df_1[,order(colSums(df_1), decreasing = T)]

df_1 <- df_1[1:50,]
vector <- rownames(df_1)

df <- df[rownames(df) %in% vector,]

df <- df[order(rowSums(df),decreasing=T),]
df <- df[,order(colSums(df), decreasing = T)]

rownames(df) <- str_remove(rownames(df), ".WP.*")
rownames(df) <- str_remove(rownames(df), "..GO.*")
rownames(df) <- str_remove(rownames(df), ".R.HSA.*")
rownames(df) <- str_replace_all(rownames(df), "\\.", " ")

rownames(df)[[48]] <- "TWEAK Signaling Pathway"

df <- as.matrix(df)
col_fun = colorRamp2(c(min(df), 0, max(df)), c("blue", "white", "red"))
col_fun(seq(min(df), max(df)))

p_scmgh_heatmap <- as.ggplot(Heatmap(df, 
                                    name = "Weights",
                                    column_title = "Proteins", row_title = "Pathways", cluster_columns = F,
                                    cluster_rows = F, col = col_fun,
                                    row_names_gp = gpar(fontsize = 7.5), 
                                    column_names_gp = gpar(fontsize = 9)
))

ggsave(sprintf("%s/Figures/Figure_5/sc1_heatmap_scmgh_ac_top20_SHAP.svg",main_dir), plot = p_scmgh_heatmap, device = "svg", width = 8, height = 7)


## FIGURE 6 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Figure_6", main_dir))

# create a super-auc with all benchmarked models 

df <- read.table(sprintf("%s/Figures/tables_auc/df_auc.csv",main_dir),
                 header = T, sep = ",")

p_super_auc <- ggplot(df, aes(x = FPR, y = TPR, fill = Model, color = Model)) +
  geom_line(aes(color=Model)) +
  labs(x = 'False Positive Rate', y = 'True Positive Rate', title = 'Receiver Operating Characteristic (ROC) Curves for Multiple Models') +
  theme_minimal()

ggsave(sprintf("%s/Figures/Figure_6/superAUC.svg",main_dir), plot = p_super_auc, device = "svg", width = 8, height = 7)

# create a lolipop plot for F1 scores

# F1 score file manually created from all models

df <- read.xlsx(sprintf("%s/Masterfiles/f1_score.xlsx", main_dir))

p_f1 <- df %>%
  arrange(F1_score) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Method=factor(Method, levels=Method)) %>%   # This trick update the factor levels
  ggplot( aes(x=Method, y=F1_score)) +
  geom_segment( aes(xend=Method, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme_bw() +
  xlab("")

ggsave(sprintf("%s/Figures/Figure_6/f1_score_lolipop.svg",main_dir), plot = p_f1, device = "svg", width = 8, height = 7)

# create heatmaps from top20 SHAP values and pathways sc1 expression PASNet output

# Mayo #

vector <- as.vector(c("IL6", "NCF2", "CCL20", "DPY30", "CHI3L1", "TNFRSF10B",
                       "POLR2F", "AREG", "CKAP4", "IL7R", "MMP8", "DCTN2", 
                       "TNFRSF6B", "LTA", "FASLG", "IL1RAP", "JUN", "ITGB1", 
                       "CD8A", "AGR2"))


df <- read.xlsx(sprintf("%s/Case_study_1/7.PASNet/Expression/Output/Mayo_test/sc1_weights_converted.xlsx",main_dir), rowNames = TRUE)
df <- df[,colnames(df) %in% vector]

df_1 <- as.data.frame(abs(df))
df_1 <- df_1[order(rowSums(df_1),decreasing=T),]
df_1 <- df_1[,order(colSums(df_1), decreasing = T)]

df_1 <- df_1[1:50,]
vector <- rownames(df_1)

df <- df[rownames(df) %in% vector,]

df <- df[order(rowSums(df),decreasing=T),]
df <- df[,order(colSums(df), decreasing = T)]

rownames(df) <- str_remove(rownames(df), ".WP.*")
rownames(df) <- str_remove(rownames(df), "..GO.*")
rownames(df) <- str_remove(rownames(df), ".R.HSA.*")
rownames(df) <- str_replace_all(rownames(df), "\\.", " ")

rownames(df)[[26]] <- "TWEAK Signaling Pathway"
rownames(df)[[28]] <- "Photodynamic therapy-induced AP-1 survival signal"

df <- as.matrix(df)
col_fun = colorRamp2(c(min(df), 0, max(df)), c("blue", "white", "red"))
col_fun(seq(min(df), max(df)))

p_mayo_heatmap_ex <- as.ggplot(Heatmap(df, 
                                    name = "Weights",
                                    column_title = "Proteins", row_title = "Pathways", cluster_columns = F,
                                    cluster_rows = F, col = col_fun,
                                    row_names_gp = gpar(fontsize = 7.5), 
                                    column_names_gp = gpar(fontsize = 9)
))

ggsave(sprintf("%s/Figures/Figure_6/sc1_heatmap_mayo_ex_top20_SHAP.svg",main_dir), plot = p_mayo_heatmap_ex, device = "svg", width = 8, height = 7)

# Stanford #

vector <- as.vector(c("CCL20", "CHI3L1", "CXCL8", "TNFRSF10B", "AREG",
                           "MMP8", "HSPA1A", "TNFRSF1A", "IL6", "CKAP4", 
                           "NCF2", "MAPK9", "LGALS9", "DPY30", "TNFRSF10A",
                           "POLR2F", "CASP1", "CLEC5A", "TNFSF10", "JUN"))

df <- read.xlsx(sprintf("%s/Case_study_1/7.PASNet/Expression/Output/Stanford_test/sc1_weights_converted.xlsx",main_dir), rowNames = TRUE)
df <- df[,colnames(df) %in% vector]

df_1 <- as.data.frame(abs(df))
df_1 <- df_1[order(rowSums(df_1),decreasing=T),]
df_1 <- df_1[,order(colSums(df_1), decreasing = T)]

df_1 <- df_1[1:50,]
vector <- rownames(df_1)

df <- df[rownames(df) %in% vector,]

df <- df[order(rowSums(df),decreasing=T),]
df <- df[,order(colSums(df), decreasing = T)]

rownames(df) <- str_remove(rownames(df), ".WP.*")
rownames(df) <- str_remove(rownames(df), "..GO.*")
rownames(df) <- str_remove(rownames(df), ".R.HSA.*")
rownames(df) <- str_replace_all(rownames(df), "\\.", " ")

rownames(df)[35] <- "TWEAK Signaling Pathway"
rownames(df)[25] <- "Photodynamic therapy-induced AP-1 survival signal"

df <- as.matrix(df)
col_fun = colorRamp2(c(min(df), 0, max(df)), c("blue", "white", "red"))
col_fun(seq(min(df), max(df)))

p_stanford_heatmap <- as.ggplot(Heatmap(df, 
                                        name = "Weights",
                                        column_title = "Proteins", row_title = "Pathways", cluster_columns = F,
                                        cluster_rows = F, col = col_fun,
                                        row_names_gp = gpar(fontsize = 7.5), 
                                        column_names_gp = gpar(fontsize = 9)
))

ggsave(sprintf("%s/Figures/Figure_6/sc1_heatmap_stanford_ex_top20_SHAP.svg",main_dir), plot = p_stanford_heatmap, device = "svg", width = 8, height = 7)


## FIGURE 7 ##
library(igraph)
library(rbioapi)

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Figure_7", main_dir))

masterfile <- read.csv(sprintf("%s/Masterfiles/APNet_masterfile.csv", main_dir))

# Nodes #
nodes_info <- read.csv(sprintf("%s/Case_study_1/5.Results/Activity/MGH/DA.csv", main_dir), header = T)
nodes_info_common <- nodes_info[nodes_info$ID %in% masterfile$ID,]
colnames(nodes_info_common)[[1]] <- "id"

# Stringify Nodes # 

Stringdb_nodes <- unique(rba_string_interactions_network(ids = nodes_info_common$id,
                                                         species = 9606,
                                                         required_score = 400))

Stringdb_nodes <- Stringdb_nodes %>% rename('source' = 'preferredName_A') %>% rename('target' = 'preferredName_B')
Stringdb_nodes <- subset(Stringdb_nodes, select = -c(stringId_A, stringId_B))

# graph_sjaracne <- graph_from_data_frame(d = network_sj_pval_common, vertices = nodes_info_common, directed = FALSE)

# write_graph(graph = graph_sjaracne, file = sprintf("%s/Figures/Figure_7/sjaracne_network.graphml", main_dir), format = "graphml")

# Driver Pathway Mapping #

pathways <- read.xlsx(sprintf("%s/Case_study_1/7.PASNet/Activity/Output/Mayo_test/sc1_weights_fixed.xlsx", main_dir), rowNames = TRUE)

pathways <- pathways[,colnames(pathways) %in% nodes_info_common$id]

data <- data.frame(matrix(NA, ncol = 3, nrow = 1))
colnames(data) <- c("source", "target", "Weight") 

for (i in 1:nrow(pathways)) {
  
  path <- pathways[rownames(pathways) %in% rownames(pathways)[i],]
  path <- as.data.frame(t(path))
  path$target <- colnames(path)[1]
  path$source <- rownames(path)
  path$Weight <- path[,1]
  path <- subset(path, select = -c(1))
  rownames(path) <- NULL
  data <- rbind(data, path)
  data <- na.omit(data)
}

pathways <- data
pathways <- filter(pathways, Weight > 0.5)

graph_pathways <- graph_from_data_frame(d = pathways, vertices = NULL, directed = FALSE)

write_graph(graph = graph_pathways, file = sprintf("%s/Figures/Figure_7/pathways_network.graphml", main_dir), 
            format = "graphml")

# Stringify SJARACNe network # 

Stringdb_nodes <- unique(rba_string_interactions_network(ids = nodes_info_common$id,
                                                         species = 9606,
                                                         required_score = 400))

Stringdb_nodes <- Stringdb_nodes %>% rename('source' = 'preferredName_A') %>% rename('target' = 'preferredName_B')
Stringdb_nodes <- subset(Stringdb_nodes, select = -c(stringId_A, stringId_B))

graph_stringdb <- graph_from_data_frame(d = Stringdb_nodes, vertices = nodes_info_top20, directed = FALSE)

write_graph(graph = graph_stringdb, file = sprintf("%s/Figures/Figure_7/stringdb_network.graphml", main_dir), 
            format = "graphml")



## SUPPLEMENTARY FIGURE 1 ##
main_dir <- "C:/Users/vasileioubill95/Desktop/Pipeline"
mkdir(sprintf("%s/Figures/Sup_Figure_1", main_dir))

# MGH #

mgh <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(mgh) <- mgh$subject_id

mgh <- mgh[,colnames(mgh) %in% c( "subject_id", "COVID", "Age_cat", "Acuity_0")]

mgh <- filter(mgh, COVID == "1")

mgh$Condition <- mgh$Acuity_0

mgh[["Condition"]] <- str_replace_all(mgh[["Condition"]], "1", "Severe")
mgh[["Condition"]] <- str_replace_all(mgh[["Condition"]], "2", "Severe")
mgh[["Condition"]] <- str_replace_all(mgh[["Condition"]], "3", "NonSevere")
mgh[["Condition"]] <- str_replace_all(mgh[["Condition"]], "4", "NonSevere")
mgh[["Condition"]] <- str_replace_all(mgh[["Condition"]], "5", "NonSevere")

mgh$Age_cat <- str_replace_all(mgh$Age_cat, "5", "80+")
mgh$Age_cat <- str_replace_all(mgh$Age_cat, "4", "65-79")
mgh$Age_cat <- str_replace_all(mgh$Age_cat, "3", "50-64")
mgh$Age_cat <- str_replace_all(mgh$Age_cat, "2", "36-49")
mgh$Age_cat <- str_replace_all(mgh$Age_cat, "1", "20-34")

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)
colnames(df) <- str_remove(colnames(df), "X")

mgh <- mgh[mgh$subject_id %in% colnames(df),]

mgh_who_score <- mgh %>%
  group_by(Acuity_0, Condition) %>%
  summarize(Sample_Count = n())

p_who <- ggplot(mgh_who_score, aes(x = Acuity_0, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "MGH - Number of Samples by WHO score",
       x = "Who Score",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

mgh_age <- mgh %>%
  group_by(Age_cat, Condition) %>%
  summarize(Sample_Count = n())

p_age <- ggplot(mgh_age, aes(x = Age_cat, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "MGH - Number of Samples by Age",
       x = "Age",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

p_mgh <- p_who | p_age

ggsave(sprintf("%s/Figures/Sup_Figure_X/MGH_who_age.svg",main_dir), plot = p_mgh, device = "svg", width = 10, height = 6)

# Mayo #

mayo <- read.xlsx(sprintf("%s/Case_study_1/1.Raw_files/Mayo/MC_raw_original.xlsx",main_dir))
rownames(mayo) <- mayo$Sample

mayo <- mayo[,colnames(mayo) %in% c( "Sample", "Diagnosis", "WHOscale", "Age_scale")]

mayo <- filter(mayo, Diagnosis == "COVID19")

mayo$Condition <- mayo$WHOscale

mayo$Condition <- str_replace_all(mayo$Condition, "1", "NonSevere")
mayo$Condition <- str_replace_all(mayo$Condition, "2", "NonSevere")
mayo$Condition <- str_replace_all(mayo$Condition, "3", "Severe")
mayo$Condition <- str_replace_all(mayo$Condition, "4", "Severe")
mayo$Condition <- str_replace_all(mayo$Condition, "5", "Severe")
mayo$Condition <- str_replace_all(mayo$Condition, "6", "Severe")
mayo$Condition <- str_replace_all(mayo$Condition, "7", "Severe")
mayo$Condition <- str_replace_all(mayo$Condition, "8", "Severe")

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Mayo/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

mayo <- mayo[mayo$Sample %in% colnames(df),]

mayo_who_score <- mayo %>%
  group_by(WHOscale, Condition) %>%
  summarize(Sample_Count = n())

p_who <- ggplot(mayo_who_score, aes(x = WHOscale, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "Mayo - Number of Samples by WHO score",
       x = "Who Score",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

mayo_age <- mayo %>%
  group_by(Age_scale, Condition) %>%
  summarize(Sample_Count = n())

p_age <- ggplot(mayo_age, aes(x = Age_scale, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "Mayo - Number of Samples by Age",
       x = "Age",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

p_mayo <- p_who | p_age

ggsave(sprintf("%s/Figures/Sup_Figure_X/Mayo_who_age.svg",main_dir), plot = p_mayo, device = "svg", width = 10, height = 6)

# Stanford #

stanford <- read.csv(sprintf("%s/Case_study_1/1.Raw_files/Stanford/PatientCharacteristics.csv",main_dir),
                  header = T)
rownames(stanford) <- stanford$sampleID

stanford <- filter(stanford, Severity != "Healthy")

stanford <- subset(stanford, select = c(sampleID ,Severity, Age_scale))

stanford$Condition <- stanford$Severity

stanford$Condition <- str_replace_all(stanford$Condition, "Asymptomatic", "NonSevere")
stanford$Condition <- str_replace_all(stanford$Condition, "Healthy", "NonSevere")
stanford$Condition <- str_replace_all(stanford$Condition, "Mild", "NonSevere")
stanford$Condition <- str_replace_all(stanford$Condition, "Moderate", "NonSevere")
stanford$Condition <- str_replace_all(stanford$Condition, "Severe", "Severe")

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Stanford/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

colnames(df) <- str_remove(colnames(df), "X")

stanford <- stanford[stanford$sampleID %in% colnames(df),]

stanford_who_score <- stanford%>%
  group_by(Severity, Condition) %>%
  summarize(Sample_Count = n())

p_who <- ggplot(stanford_who_score, aes(x = Severity, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "Stanford - Number of Samples by WHO score",
       x = "Who Score",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

stanford_age <- stanford %>%
  group_by(Age_scale, Condition) %>%
  summarize(Sample_Count = n())

p_age <- ggplot(stanford_age, aes(x = Age_scale, y = Sample_Count, fill = Condition)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(aes(label = Sample_Count), position = position_dodge(width = 1), vjust = -0.5) +  # Add labels
  labs(title = "Stanford - Number of Samples by Age",
       x = "Age",
       y = "Number of Samples",
       fill = "Condition") +
  theme_minimal()

p_stanford <- p_who | p_age

ggsave(sprintf("%s/Figures/Sup_Figure_X/Stanford_who_age.svg",main_dir), plot = p_stanford, device = "svg", width = 10, height = 6)



## SUPPLEMENTARY FIGURE 2 ##
library(NetBID2)

rm(list=ls())
# load directories
main_dir <- "C:/Users/vasileioubill95/Desktop/Pipeline"
mkdir(sprintf("%s/Figures/Sup_Figure_2", main_dir))

dataset <- "Stanford" # MGH / Mayo / Stanford each time

## load NetBID2 results
comp_name <- 'Severe.Vs.NonSevere'
sit_1 <- 'Severe'
sit_0 <- 'NonSevere'
Condition <- "Condition"

network.dir <- sprintf("%s/Case_study_1/3.NetBID2/%s/project", main_dir, dataset)
analysis.par <- list()
analysis.par$out.dir.DATA <- sprintf("%s/Case_study_1/3.NetBID2/%s/Driver_output/driver_project/DATA", 
                                     main_dir, 
                                     dataset) 

NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

DA <- analysis.par$DA$Severe.Vs.NonSevere

DA$ID <- str_replace_all(DA$ID, "_SIG", "")
DA$ID <- str_replace_all(DA$ID, "_TF", "")

DA <- DA[!duplicated(DA[,c("ID")]),]
rownames(DA) <- DA$ID

DA_up <- filter(DA, logFC > 0)
driver_list_up <- (DA_up[order(DA_up$adj.P.Val,decreasing=FALSE),])[1:5,]
DA_down <- filter(DA , logFC < 0)
driver_list_down <- (DA_down[order(DA_down$adj.P.Val,decreasing=FALSE),])[1:5,]

driver_list <- rbind(driver_list_up, driver_list_down)
driver_list <- driver_list$ID

DE <- analysis.par$DE$Severe.Vs.NonSevere
target_list <- analysis.par$tf.network$target_list
names(target_list) <- str_remove(names(target_list), "_TF")

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=driver_list,
                 driver_DA_Z=DA[driver_list,'Z-statistics'],
                 driver_DE_Z=DE[driver_list,'Z-statistics'],
                 target_list=target_list,
                 top_driver_number=10,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in Severe',right_annotation = 'high in NonSevere',
                 main=sprintf("%s Top 10 Drivers", dataset),target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf("%s/Figures/Sup_Figure_2/%s_sig_drivers.pdf", main_dir, dataset))

gs.preload(use_spe='Homo sapiens',update=FALSE)
res <- funcEnrich.Fisher(input_list=driver_list,bg_list=unique(DA$ID),use_gs=c('H','CP:REACTOME','BP','CGP'), Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)

draw.funcEnrich.cluster(funcEnrich_res= res,top_number=10,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf("%s/Figures/Sup_Figure_2/%s_sig_drivers_heatmap.pdf", main_dir, dataset),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)

## SUPPLEMENTARY FIGURE 3 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Sup_Figure_3", main_dir))

# create bubble plot for DEPs and Hiddens proteins by using REACTOME db

DEPs_paths <- read.table(sprintf("%s/Case_study_1/EnrichR_KG/DEPs_Hidden/Reactome_2022_common_DEPs.txt", main_dir), header=T, sep = "\t")

Hidden_paths <- read.table(sprintf("%s/Case_study_1/EnrichR_KG/DEPs_Hidden/Reactome_2022_common_Hidden.txt", main_dir), header=T, sep = "\t")

DEPs_paths <- filter(DEPs_paths, Adjusted.P.value < 0.05)
DEPs_paths <- DEPs_paths[order(DEPs_paths$Adjusted.P.value, decreasing = FALSE), ]

DEPs_paths <- DEPs_paths[1:30,]

DEPs_paths$logqvalue <- -log10(DEPs_paths$Adjusted.P.value)
DEPs_paths <- DEPs_paths[order(DEPs_paths$logqvalue, decreasing = F),]

DEPs_paths$Term <- str_remove(DEPs_paths$Term, ".WP.*")
DEPs_paths$Term <- str_remove(DEPs_paths$Term, "..GO.*")
DEPs_paths$Term <- str_remove(DEPs_paths$Term, ".R.HSA.*")
DEPs_paths$Term <- str_replace_all(DEPs_paths$Term, "\\.", " ")

DEPs_paths$Term[10] <- "Diseases Of Signal Transduction By GFR And Second Messengers"
DEPs_paths$Term[18] <- "Immunoregulation Between A Lymphoid And A non-Lymphoid Cell"

# Plot
p_bubble_deps <- ggplot(data = DEPs_paths, aes(x=logqvalue, y=Term, size = Combined.Score , color = Odds.Ratio )) +
  geom_point(alpha=1) +
  ggtitle("Reactome - DEPs") +
  ylab("Terms") +
  xlab("-10log(P Value)") +
  labs(color = "Z Score") +
  labs(size = "Combine Score") +
  scale_y_discrete(limits = DEPs_paths$Term) +
  scale_colour_continuous(trans = 'reverse') +
  scale_color_gradient(low = "firebrick1", high = "firebrick")

ggsave(sprintf("%s/Figures/Sup_Figure_1/bubble_plot_deps.svg",main_dir), plot = p_bubble_deps, device = "svg", width = 10, height = 7)

Hidden_paths <- filter(Hidden_paths, Adjusted.P.value < 0.05)
Hidden_paths <- Hidden_paths[order(Hidden_paths$Adjusted.P.value, decreasing = FALSE), ]

Hidden_paths <- Hidden_paths[1:30,]

Hidden_paths$logqvalue <- -log10(Hidden_paths$Adjusted.P.value)
Hidden_paths <- Hidden_paths[order(Hidden_paths$logqvalue, decreasing = F),]

Hidden_paths$Term <- str_remove(Hidden_paths$Term, ".WP.*")
Hidden_paths$Term <- str_remove(Hidden_paths$Term, "..GO.*")
Hidden_paths$Term <- str_remove(Hidden_paths$Term, ".R.HSA.*")
Hidden_paths$Term <- str_replace_all(Hidden_paths$Term, "\\.", " ")

Hidden_paths$Term[11] <- "Immunoregulation Between A Lymphoid And A non-Lymphoid Cell"

# Plot
p_bubble_hidden <- ggplot(data = Hidden_paths, aes(x=logqvalue, y=Term, size = Combined.Score , color = Odds.Ratio )) +
  geom_point(alpha=1) +
  ggtitle("Reactome - Hidden") +
  ylab("Terms") +
  xlab("-10log(P Value)") +
  labs(color = "Z Score") +
  labs(size = "Combine Score") +
  scale_y_discrete(limits = Hidden_paths$Term) +
  scale_colour_continuous(trans = 'reverse') +
  scale_color_gradient(low = "royalblue1", high = "royalblue4")

ggsave(sprintf("%s/Figures/Sup_Figure_1/bubble_plot_hidden.svg",main_dir), plot = p_bubble_hidden, device = "svg", width = 10, height = 7)


## SUPPLEMENTARY FIGURE 4 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Sup_Figure_4", main_dir))

# GEORGE CODE

## SUPPLEMENTARY FIGURE 5 ##

# create VennDiagrams from common files

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Sup_Figure_3", main_dir))

# create VennDiagrams from common files

mgh_neg <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/MGH/MGH_neg.txt", main_dir), header = T)
mgh_neg <- mgh_neg$x
mgh_pos <- read.table(sprintf("%s/Case_study_1/6.Overlaps/Activity/MGH/MGH_pos.txt", main_dir), header = T)
mgh_pos <- mgh_pos$x

scmgh_neg <- read.table(sprintf("%s/Case_study_2/5.Overlaps/scMGH_neg.txt", main_dir), header = T)
scmgh_neg <- scmgh_neg$x
scmgh_pos <- read.table(sprintf("%s/Case_study_2/5.Overlaps/scMGH_pos.txt", main_dir), header = T)
scmgh_pos <- scmgh_pos$x

pos <- list(MGH=sample(mgh_pos),
            scMGH=sample(scmgh_pos))

neg <- list(MGH=sample(mgh_neg),
            scMGH=sample(scmgh_neg))

p_pos <- ggvenn::ggvenn(pos, set_name_size = 4, text_size = 4, show_percentage = FALSE) + 
  labs(title = "Common Positive DA") +
  theme(text = element_text(size=10))

p_neg <-  ggvenn::ggvenn(neg, set_name_size = 4, text_size = 4, show_percentage = FALSE) + 
  labs(title = "Common Negative DA") +
  theme(text = element_text(size=10))

p_venn <- p_pos + p_neg

ggsave(sprintf("%s/Figures/Sup_Figure_3/venndiagram.svg",main_dir), plot = p_venn, device = "svg", width = 8, height = 7)

enrich_pos <- read.csv(sprintf("%s/Case_study_2/EnrichR_KG/Activity/Enrichr-KG_pos.csv",main_dir))
enrich_pos <- enrich_pos[order(enrich_pos$q.value, decreasing = FALSE), ]

enrich_pos$Term <- str_remove(enrich_pos$Term, ".WP.*")
enrich_pos$Term <- str_remove(enrich_pos$Term, "..GO.*")
enrich_pos$Term <- str_remove(enrich_pos$Term, ".R.HSA.*")
enrich_pos$Term <- str_replace_all(enrich_pos$Term, "\\.", " ")

enrich_pos <- enrich_pos[!duplicated(enrich_pos[,c("Term")]),]

enrich_pos <- enrich_pos[1:30,]

enrich_pos$logqvalue <- -log10(enrich_pos$q.value)
enrich_pos <- enrich_pos[order(enrich_pos$logqvalue, decreasing = F),]

## SUPPLEMENTARY FIGURE 6 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Sup_Figure_6", main_dir))

# create heatmap with top20 SHAP values from PASNet by using MGH expression matrix

# Mayo SHAP #

vector <- as.vector(c("IL6", "NCF2", "CCL20", "DPY30", "CHI3L1", "TNFRSF10B",
                       "POLR2F", "AREG", "CKAP4", "IL7R", "MMP8", "DCTN2", 
                       "TNFRSF6B", "LTA", "FASLG", "IL1RAP", "JUN", "ITGB1", 
                       "CD8A", "AGR2"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Expression/MGH/count_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH PASNet - Mayo top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_6/mgh_mayo_heatmap_expression.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford SHAP #

vector <- as.vector(c("CCL20", "CHI3L1", "CXCL8", "TNFRSF10B", "AREG",
                   "MMP8", "HSPA1A", "TNFRSF1A", "IL6", "CKAP4", 
                   "NCF2", "MAPK9", "LGALS9", "DPY30", "TNFRSF10A",
                   "POLR2F", "CASP1", "CLEC5A", "TNFSF10", "JUN"))


df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Expression/MGH/count_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH PASNet - Stanford top20 SHAP")

ggsave(sprintf("%s/Figures/Sup_Figure_6/mgh_stanford_heatmap_expression.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# create heatmap with top20 SHAP values from PASNet by using Mayo and Stanford expression matrices

# Mayo #

vector <- as.vector(c("IL6", "NCF2", "CCL20", "DPY30", "CHI3L1", "TNFRSF10B",
                      "POLR2F", "AREG", "CKAP4", "IL7R", "MMP8", "DCTN2", 
                      "TNFRSF6B", "LTA", "FASLG", "IL1RAP", "JUN", "ITGB1", 
                      "CD8A", "AGR2"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Expression/Mayo/count_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]

pheno <- read.xlsx(sprintf("%s/Case_study_1/1.Raw_files/Mayo/MC_raw.xlsx",main_dir))
rownames(pheno) <- pheno$Sample

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]

df <- df[,colnames(df) %in% rownames(pheno)]

matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- subset(pheno, select = c(WHOscale))
colnames(pheno)[colnames(pheno) == "WHOscale"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "1", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "2", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "3", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "4", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "5", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "6", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "7", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "8", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(Severe="red", NonSevere="blue"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "Mayo PASNet - Mayo top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_6/mayo_heatmap_expression.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford #

vector <- as.vector(c("CCL20", "CHI3L1", "CXCL8", "TNFRSF10B", "AREG",
                      "MMP8", "HSPA1A", "TNFRSF1A", "IL6", "CKAP4", 
                      "NCF2", "MAPK9", "LGALS9", "DPY30", "TNFRSF10A",
                      "POLR2F", "CASP1", "CLEC5A", "TNFSF10", "JUN"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Expression/Stanford/count_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")

pheno <- read.csv(sprintf("%s/Case_study_1/1.Raw_files/Stanford/PatientCharacteristics.csv",main_dir),
                  header = T)
rownames(pheno) <- pheno$sampleID

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]
pheno <- pheno[order(rownames(pheno)),]

matrix <- as.matrix(df)

matrix_norm <- scale(t(matrix), center = T, scale = T)
matrix_norm <- matrix_norm[order(rownames(matrix_norm)),]
matrix_norm <- t(matrix_norm)

pheno <- subset(pheno, select = c(Severity))
colnames(pheno)[colnames(pheno) == "Severity"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "Asymptomatic", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Healthy", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Mild", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Moderate", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Severe", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(Severe="red", NonSevere="blue"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "Stanford PASNet - Stanford top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_6/stanford_heatmap_expression.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

## SUPPLEMENTARY FIGURE 7 ##

rm(list=setdiff(ls(), "main_dir"))
gc()
mkdir(sprintf("%s/Figures/Sup_Figure_7", main_dir))

# create heatmap with top20 SHAP values from Random Forest by using MGH activity matrix

# Mayo SHAP #

vector <- as.vector(c("POLR2F", "CCL22", "DDAH1", "FKBP5", "MSR1", 
                      "SERPINB8", "NCF2", "KRT18", "MAPK9", "MASP1", 
                      "LBR", "FMNL1", "GPNMB", "CST7", "KLK8", "GFRA2", 
                      "CD207", "ICOSLG", "ASAH2", "F7"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH Activity RF - Mayo top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_7/mgh_mayo_heatmap_activity_rf.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford SHAP #

vector <- as.vector(c("POLR2F", "CCL22", "DDAH1", "FKBP5", "MSR1", 
                      "SERPINB8", "MAPK9", "KRT18", "MASP1", "FMNL1",
                      "NCF2", "LBR", "GPNMB", "CST7", "ICOSLG", "CD207",
                      "KLK8", "GFRA2", "TNFSF10", "ASAH2"))


df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH Activity RF - Stanford top20 SHAP")

ggsave(sprintf("%s/Figures/Sup_Figure_7/mgh_stanford_heatmap_activity_rf.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# scMGH SHAP #

vector <- as.vector(c("KRT18", "F3", "FLT3LG", "BAG3","CD48", "ATP6AP2", "S100A11", "SLAMF6", 
                 "MAECO", "JUN", "VCAN", "IL5RA", "CD1C", "NCF2", "SPP1", "CD4",
                 "ENAH", "MDK", "FEN1", "MAPK9"))


df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/MGH/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")
matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- read.table(sprintf("%s/Case_study_1/1.Raw_files/MGH/MGH_COVID_Clinical_Info.txt",main_dir), header = T, sep = ";")

rownames(pheno) <- pheno$subject_id
pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[pheno$subject_id %in% colnames(matrix_norm),]
pheno <- pheno[,colnames(pheno) %in% c("KIDNEY", "DIABETES", "Age_cat", "HEART", "Acuity_0")]

colnames(pheno)[colnames(pheno) == "Acuity_0"] = "Condition"

pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "1", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "2", "Severe")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "3", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "4", "NonSevere")
pheno[["Condition"]] <- str_replace_all(pheno[["Condition"]], "5", "NonSevere")

pheno$HEART <- str_replace_all(pheno$HEART, "0", "no")
pheno$HEART <- str_replace_all(pheno$HEART, "1", "yes")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "0", "no")
pheno$KIDNEY <- str_replace_all(pheno$KIDNEY, "1", "yes")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "0", "no")
pheno$DIABETES <- str_replace_all(pheno$DIABETES, "1", "yes")

pheno$Age_cat <- str_replace_all(pheno$Age_cat, "5", "80+")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "4", "65-79")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "3", "50-64")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "2", "36-49")
pheno$Age_cat <- str_replace_all(pheno$Age_cat, "1", "20-34")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        Diabetes = pheno$DIABETES,
                        Heart = pheno$HEART, 
                        Kidney = pheno$KIDNEY,
                        Age = pheno$Age_cat,
                        col = list(Condition = c(Severe="red", NonSevere="blue"), 
                                   Heart = c("no" = "purple", "yes" = "darkgrey"),
                                   Kidney = c("no" = "yellow", "yes" = "deeppink"),
                                   Diabetes = c("no"= "green", "yes" = "brown"),
                                   Age = c("20-34" = "lightgreen", "36-49" = "yellowgreen","50-64"= "darkgreen", "65-79"="blue2","80+"= "blue4"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "MGH Activity RF - scMGH top20 SHAP")

ggsave(sprintf("%s/Figures/Sup_Figure_7/mgh_scmgh_heatmap_activity_rf.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)


# create heatmap with top20 SHAP values from Random Forest by using Mayo and Stanford activity matrices

# Mayo #

vector <- as.vector(c("POLR2F", "CCL22", "DDAH1", "FKBP5", "MSR1", 
                      "SERPINB8", "NCF2", "KRT18", "MAPK9", "MASP1", 
                      "LBR", "FMNL1", "GPNMB", "CST7", "KLK8", "GFRA2", 
                      "CD207", "ICOSLG", "ASAH2", "F7"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Mayo/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]

pheno <- read.xlsx(sprintf("%s/Case_study_1/1.Raw_files/Mayo/MC_raw.xlsx",main_dir))
rownames(pheno) <- pheno$Sample

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]

df <- df[,colnames(df) %in% rownames(pheno)]

matrix <- as.matrix(df)

matrix_norm <- t(scale(t(matrix), center = T, scale = T))

pheno <- subset(pheno, select = c(WHOscale))
colnames(pheno)[colnames(pheno) == "WHOscale"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "1", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "2", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "3", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "4", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "5", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "6", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "7", "Severe")
pheno$Condition <- str_replace_all(pheno$Condition, "8", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(Severe="red", NonSevere="blue"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "Mayo Activity RF - Mayo top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_7/mayo_heatmap_activity_rf.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

# Stanford #

vector <- as.vector(c("POLR2F", "CCL22", "DDAH1", "FKBP5", "MSR1", 
                      "SERPINB8", "MAPK9", "KRT18", "MASP1", "FMNL1",
                      "NCF2", "LBR", "GPNMB", "CST7", "ICOSLG", "CD207",
                      "KLK8", "GFRA2", "TNFSF10", "ASAH2"))

df <- read.csv(sprintf("%s/Case_study_1/4.Matrices/Activity/Stanford/activity_matrix.csv",main_dir),
               header = T,
               row.names = 1)

df <- df[rownames(df) %in% vector,]
colnames(df) <- str_remove(colnames(df), "X")

pheno <- read.csv(sprintf("%s/Case_study_1/1.Raw_files/Stanford/PatientCharacteristics.csv",main_dir),
                  header = T)
rownames(pheno) <- pheno$sampleID

pheno <- pheno %>% mutate(across(everything(), as.character))
pheno <- pheno[rownames(pheno) %in% colnames(df),]
pheno <- pheno[order(rownames(pheno)),]

matrix <- as.matrix(df)

matrix_norm <- scale(t(matrix), center = T, scale = T)
matrix_norm <- matrix_norm[order(rownames(matrix_norm)),]
matrix_norm <- t(matrix_norm)

pheno <- subset(pheno, select = c(Severity))
colnames(pheno)[colnames(pheno) == "Severity"] = "Condition"

pheno$Condition <- str_replace_all(pheno$Condition, "Asymptomatic", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Healthy", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Mild", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Moderate", "NonSevere")
pheno$Condition <- str_replace_all(pheno$Condition, "Severe", "Severe")

ha <- HeatmapAnnotation(Condition = pheno$Condition,
                        col = list(Condition = c(Severe="red", NonSevere="blue"))
)

p <- Heatmap(matrix_norm, km = 1, 
             show_row_names = T, 
             show_column_names = F, 
             cluster_rows = T, 
             cluster_columns = T, 
             name = "Activity", 
             top_annotation = ha,
             column_title = "Stanford Activity RF - Stanford top20 SHAP")


ggsave(sprintf("%s/Figures/Sup_Figure_7/stanford_heatmap_activity_rf.svg",main_dir), plot = as.ggplot(p), device = "svg", width = 6, height = 5)

