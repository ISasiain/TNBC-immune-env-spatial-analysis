#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(tidyr)
library(ComplexHeatmap)
library(survivalAnalysis)
library(survival)
library(survminer)

mIHC_counts1 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_1_dataframe.csv")
mIHC_counts2 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_2_dataframe.csv")
mIHC_counts3 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_3_dataframe.csv")
mIHC_counts4 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_4_dataframe.csv")
mIHC_counts5 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_5_dataframe.csv")


mIHC_counts <- rbind(mIHC_counts1, mIHC_counts2, mIHC_counts3, mIHC_counts4, mIHC_counts5)


annotation_file <- read.csv2("/home/isc/Spatial_immune_env/data_from_suze/data/supplData_withimages.csv")
rownames(annotation_file) <- annotation_file$uid

#
# COMPARING CELL COUNTS WITH RNAseq DATA
#


# Comparing CD68 counts with RNAseq data

vector_of_counts <- as.numeric(mIHC_counts$CD68)
vector_of_mRNAseq <- as.numeric(sapply(mIHC_counts$TMArQ_CORE_ID, function(core) annotation_file[core, "CibersortX.macrophage"]))

ggplot(data=data.frame(vector_of_counts, vector_of_mRNAseq), aes(x=log(vector_of_counts), y=vector_of_mRNAseq)) +
  geom_point() +
  xlab("log Macrophague counts") +
  ylab("Cibersort Macrophague fraction") +
  theme_classic()

# Comparing PAN-CK counts with RNAseq data

vector_of_counts <- as.numeric(mIHC_counts$PAN.CK)
vector_of_mRNAseq <- as.numeric(sapply(mIHC_counts$TMArQ_CORE_ID, function(core )annotation_file[core, "CibersortX.epithelial"]))
vector_of_ASCAT <- as.numeric(sapply(mIHC_counts$TMArQ_CORE_ID, function(core )annotation_file[core, "ASCAT_TUM_FRAC"]))


ggplot(data=data.frame(vector_of_counts, vector_of_mRNAseq), aes(x=log(vector_of_counts), y=vector_of_mRNAseq)) +
  geom_point() +
  xlab("log Tumour cell counts") +
  ylab("Cibersort Epithelial fraction") +
  theme_classic()


ggplot(data=data.frame(vector_of_counts, vector_of_ASCAT), aes(x=log(vector_of_counts), y=vector_of_ASCAT)) +
  geom_point() +
  xlab("log Tumour cell counts") +
  ylab("ASCAT tumour purity") +
  theme_classic()


# Comparing CD20 counts with RNAseq data

vector_of_counts <- as.numeric(mIHC_counts$CD20)
vector_of_mRNAseq <-as.numeric(sapply(mIHC_counts$TMArQ_CORE_ID, function(core )annotation_file[core, "CibersortX.Bcell"]))

ggplot(data=data.frame(vector_of_counts, vector_of_mRNAseq), aes(x=log(vector_of_counts), y=vector_of_mRNAseq)) +
  geom_point() +
  xlab("log CD20 cell counts") +
  ylab("Cibersort B lymophocyte fraction") +
  theme_classic()

#
# CD20 VS PATHOLOGIST
#

# Vector of CD20 paths
vector_of_counts <- log(mIHC_counts$CD20 +1) # Adding one to avoid getting -Inf
vector_of_CD20_path <- as.factor(unname(sapply(mIHC_counts$TMArQ_CORE_ID, function(core) annotation_file[core, "CD20path"])))

#General kruskal willis test
kruskal_res <- kruskal.test(vector_of_counts ~ vector_of_CD20_path)

#Wilcox test

pairwise_res <- compare_means(data = data.frame(vector_of_counts, vector_of_CD20_path),
                              vector_of_counts ~ vector_of_CD20_path,
                              method = "wilcox.test")
my_comparison <- list(c("0", "1"),c("1","2"), c("2","3"))


ggplot(data=na.omit(data.frame(vector_of_counts, vector_of_CD20_path)), aes(x=vector_of_CD20_path, y=vector_of_counts)) +
  geom_violin(width=1, aes(fill=vector_of_CD20_path)) +
  geom_jitter(size=0.75, width=0.08) +
  geom_boxplot(width=0.075) +
  labs(fill="Pathologist's score", x = "CD20 pathologist's score", y = "log B lymphocyte counts") +
  annotate("text", x=3.75, y=0.6, label=paste0("Kruskal-Wallis p value:\n", signif(kruskal_res$p.value, 2)), size = 7) + # Adjust size here
  stat_compare_means(pairwise_res, comparisons = my_comparison, size=7) +
  theme_classic()

        


# Analysing the effect of cell counts in survival

# Assuming annotation_file and mIHC_counts are your data frames

# Subset data
time <- annotation_file[mIHC_counts$TMArQ_CORE_ID, "time"]
status <- annotation_file[mIHC_counts$TMArQ_CORE_ID, "event"]
tils <- annotation_file[mIHC_counts$TMArQ_CORE_ID, "TILs"]

# Ensure 'time' and 'status' are of the same length and type
# Check if any NA values are present and handle them if necessary

# Combine into a data frame
surv_data <- data.frame(time = as.numeric(time), status = status,
                         TILs = as.numeric(tils),
                         PAN.CK = log(mIHC_counts$PAN.CK  +1),
                         CD8 = log(mIHC_counts$CD8 +1),
                         CD4 = log(mIHC_counts$CD4 +1),
                         CD20 = log(mIHC_counts$CD20 +1),
                         CD68 = log(mIHC_counts$CD68 +1),
                         CD4_FOXP3 = log(mIHC_counts$CD4_FOXP3 +1),
                         CD8_FOXP3 = log(mIHC_counts$CD8_FOXP3 +1),
                         Other = log(mIHC_counts$Other +1))

# Perform survival analysis
result = survivalAnalysis::analyse_multivariate(data = surv_data,
                     time_status = c("time", "status"),
                     covariates = c("TILs","CD8", "CD4", "CD20", "CD68", "CD4_FOXP3", "CD8_FOXP3"))

survivalAnalysis::forest_plot(result)


#Plotting survival plot of Treg infiltration in high TILs cores

surv_data2 <- surv_data[!is.na(surv_data$TILs) & surv_data$TILs > 30,]                      

# Effect of T cyt
cutoff <- median(log(mIHC_counts$CD8+1))
surv_data2$Tcyt_status <- as.factor(sapply(surv_data2$CD8, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tcyt_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tcyt_status <- as.factor(sapply(surv_data2$CD8, function(val) {if (val > 6.2) {"high"} else if (val > 4.45) {"medium"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tcyt_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tcyt_status <- as.factor(sapply(surv_data2$CD8, function(val) {if (val > 6.7) {"high"} else if (val > 5.5) {"medium-high"} else if (val > 3.8) {"medium-low"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tcyt_status, data = surv_data2), risk.table = TRUE)


# Effect of T helper

cutoff <- median(log(mIHC_counts$CD4))
surv_data2$Th_status <- as.factor(sapply(surv_data2$CD4, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Th_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Th_status <- as.factor(sapply(surv_data2$CD4, function(val) {if (val > 6.25) {"high"} else if (val > 4.05) {"medium"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Th_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Th_status <- as.factor(sapply(surv_data2$CD4, function(val) {if (val > 6.6) {"high"} else if (val > 5.5) {"medium-high"} else if (val > 3.8) {"medium-low"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Th_status, data = surv_data2), risk.table = TRUE)

# Effect of cd4+ T regulatory

cutoff <- median(log(mIHC_counts$CD4_FOXP3))
surv_data2$Tregcd4_status <- as.factor(sapply(surv_data2$CD4_FOXP3, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd4_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd4_status <- as.factor(sapply(surv_data2$CD4_FOXP3, function(val) {if (val > 4.8) {"high"} else if (val > 3) {"medium"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd4_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd4_status <- as.factor(sapply(surv_data2$CD4_FOXP3, function(val) {if (val > 5.1) {"high"} else if (val > 4) {"medium-high"} else if (val > 2.3) {"medium-low"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd4_status, data = surv_data2), risk.table = TRUE)


# Effect of cd8+ T regulatory

cutoff <- median(log(mIHC_counts$CD8_FOXP3))
surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > 3) {"high"} else if (val > 0.6) {"medium"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > 3.9) {"high"} else if (val > 2.1) {"medium-high"} else if (val > 0.1) {"medium-low"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE)


# Effect of other cells
cutoff <- median(log(mIHC_counts$Other))
surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > 1.5) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > 3) {"high"} else if (val > 0.6) {"medium"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE, conf.int = TRUE)


surv_data2$Tregcd8_status <- as.factor(sapply(surv_data2$CD8_FOXP3, function(val) {if (val > 3.9) {"high"} else if (val > 2.1) {"medium-high"} else if (val > 0.1) {"medium-low"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tregcd8_status, data = surv_data2), risk.table = TRUE)





ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Treg_status, data = surv_data2))

summary(surv_data2$Treg_status)

plot(log(mIHC_counts$CD8), log(mIHC_counts$CD))

# Plotting barplot with cell counts

mIHC_counts$TILs <- as.numeric(
  
  annotation_file[mIHC_counts$TMArQ_CORE_ID, "TILs"]
  
)


# Convert CORE_ID to factor for correct ordering
mIHC_counts$CORE_ID <- factor(mIHC_counts$CORE_ID, levels = mIHC_counts$CORE_ID[order(as.numeric(annotation_file[mIHC_counts$TMArQ_CORE_ID, "TILs"]))])


# Adding cluster annotation to barplot and hearmap

dins <- readRDS("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/DIN_matrices/r100_all_markers.rds")
cls <- readRDS("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/detected_clusters/r100_allclusters_clustered_cells.rds")
anno <- readRDS("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/detected_clusters/r100_allclusters.rds")


#Generating lsit to store 
cls_ls <- list()

# Determine which samples have at least one CD20 immune cluster
for (marker in c("CD20", "CD8", "CD4")) {
  
  for (core in names(anno)) {
    
    if (!"Low cell density" %in% anno[[core]]) {
      
      # Getting clusters annotated as tumour
      if (anno[[core]][["type"]] == "Heterogeneous") {
        
        for (cluster in 1:length(anno[[core]]$clusters)) {
          
          if (anno[[core]][["clusters"]][[cluster]][["type"]] == "IMMUNE") {
            
            if (names(which.max(table(dins[[core]][cls[[core]]$cluster == cluster, "Phenotype"])[c("CD20", "CD4", "CD8")])) == marker) {
              
              # Append the core value to the list
              cls_ls[[paste0(marker)]] <- c(cls_ls[[paste0(marker)]], core)
              
            }
            
          } else if (anno[[core]][["clusters"]][[cluster]][["type"]] == "TUMOUR") {
            
            cls_ls[["Tumour"]] <- c(cls_ls[["Tumour"]], core)
             
          } else if (anno[[core]][["clusters"]][[cluster]][["type"]] == "MIXED") {
            
            cls_ls[["Mixed"]] <- c(cls_ls[["Mixed"]], core)
            
          }
        }
      } else if (anno[[core]][["type"]] == "Homogeneous") {
        
        if (anno[[core]][["class"]][["type"]] == "TUMOUR") {
          
          cls_ls[["Homo_tum"]] <- c(cls_ls[["Homo_tum"]], core)
        
        } else if (anno[[core]][["class"]][["type"]] == "IMMUNE") {
        
          cls_ls[["Homo_immune"]] <- c(cls_ls[["Homo_immune"]], core)
          
        }  else if (anno[[core]][["class"]][["type"]] == "MIXED") {
          
          cls_ls[["Homo_mixed"]] <- c(cls_ls[["Homo_mixed"]], core)
          
        }
    }
  }
  }
}


vec_of_blocks_cd4 <- lapply(cls_ls$CD4, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_cd8 <- lapply(cls_ls$CD8, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_cd20 <- lapply(cls_ls$CD20, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_tum <- lapply(cls_ls$Tumour, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_mix <- lapply(cls_ls$Mixed, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_htum <- lapply(cls_ls$Homo_tum, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_himm <- lapply(cls_ls$Homo_immune, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))

vec_of_blocks_hmix <- lapply(cls_ls$Homo_mixed, function(path) sub(".*/Core_(\\d+)_(\\w+)(\\d+)_spe\\.object\\.rds", "BLOCK_\\1|\\2\\3", path))


# Adding information to mIHC_counts
mIHC_counts$CD20_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_cd20) {1} else {0}})
mIHC_counts$CD8_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_cd8) {1} else {0}})
mIHC_counts$CD4_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_cd4) {1} else {0}})
mIHC_counts$Tum_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_tum) {1} else {0}})
mIHC_counts$Mix_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_mix) {1} else {0}})
mIHC_counts$HTum_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_htum) {1} else {0}})
mIHC_counts$HImm_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_himm) {1} else {0}})
mIHC_counts$HMix_clusters <- sapply(mIHC_counts$CORE_ID, function(id) {if (id %in% vec_of_blocks_hmix) {1} else {0}})




# Gather the dataframe
dataframe_long <- mIHC_counts %>% 
  gather(key = "Phenotype", value = "Count", c(6:12,))

# Convert Count to numeric
dataframe_long$Count <- as.numeric(dataframe_long$Count)

# Create the plot
p1 <- ggplot(data = dataframe_long[!is.na(dataframe_long$TILs),], aes(x = CORE_ID, fill = Phenotype, y = Count)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = -25)) +
  xlab("Cores") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 0))

# Create the second plot with filled bars
p2 <- ggplot(data = dataframe_long[!is.na(dataframe_long$TILs),], aes(x = CORE_ID, fill = Phenotype, y = Count)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = -25)) +
  xlab("Cores") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 0))

# Print the plots
print(p1)

print(p2)



# Adding TILs to barplot

p1 + geom_point(aes(x=CORE_ID, y=TILs*100)) +
  scale_y_continuous(name = "Cell counts",sec.axis = sec_axis(name = "TILs",  ~.*1/100))


p2 + geom_point(aes(x=CORE_ID, y=1-TILs/100)) +
     scale_y_continuous(name = "Cell counts",sec.axis = sec_axis(name = "TILs",  ~.*100 - 100))


# Adding annotation. Samples with CD20 immune clusters
p1 + annotate("segment", x=mIHC_counts[mIHC_counts$CD20_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD20_clusters == 1, "CORE_ID"], y =-500, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
p2 + annotate("segment", x=mIHC_counts[mIHC_counts$CD20_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD20_clusters == 1, "CORE_ID"], y =-0.04, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))

# Adding annotation. Samples with CD8 immune clusters
p1 + annotate("segment", x=mIHC_counts[mIHC_counts$CD8_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD8_clusters == 1, "CORE_ID"], y =-500, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
p2 + annotate("segment", x=mIHC_counts[mIHC_counts$CD8_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD8_clusters == 1, "CORE_ID"], y =-0.04, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))

# Adding annotation. Cores with CD4 immune clusters
p1 + annotate("segment", x=mIHC_counts[mIHC_counts$CD4_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD4_clusters == 1, "CORE_ID"], y =-500, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
p2 + annotate("segment", x=mIHC_counts[mIHC_counts$CD4_clusters == 1, "CORE_ID"], xend = mIHC_counts[mIHC_counts$CD4_clusters == 1, "CORE_ID"], y =-0.04, yend = 0, arrow = arrow(length = unit(0.06, "inches"), type = "closed"))


# Generating heatmap

mIHC_counts$TILs <- as.numeric(annotation_file[mIHC_counts$TMArQ_CORE_ID, "TILs"])
mIHC_counts$IMstatus <- as.factor(annotation_file[mIHC_counts$TMArQ_CORE_ID, "TNBCtype_IMpositive"])
mIHC_counts$ASCAT <- as.numeric(annotation_file[mIHC_counts$TMArQ_CORE_ID, "ASCAT_TUM_FRAC"])
mIHC_counts$HRD <- as.factor(annotation_file[mIHC_counts$TMArQ_CORE_ID, "HRD.2.status"])
mIHC_counts$PAM50 <- as.factor(annotation_file[mIHC_counts$TMArQ_CORE_ID, "PAM50_NCN"])
mIHC_counts$subtypes_7 <- as.factor(annotation_file[mIHC_counts$TMArQ_CORE_ID, "TNBCtype_org"])
mIHC_counts$subtypes_4 <- as.factor(annotation_file[mIHC_counts$TMArQ_CORE_ID, "TNBCtype4_n235_notPreCentered"])

# Merging PAM50 annotations
mIHC_counts$PAM50 <- sapply(mIHC_counts$PAM50, function(val) {if (is.na(val)) {NA}
                                                              else if (val=="Basal") {"Basal"} 
                                                              else if (val == "unclassified") {"Unclassified"}
                                                              else {"Non Basal"}})


annotations_left <- rowAnnotation(TILs = mIHC_counts$TILs,
                                  HRD_status = mIHC_counts$HRD,
                                  PAM50_subtype = mIHC_counts$PAM50,
                                  IM_status = mIHC_counts$IMstatus,
                                  Refined_subtypes = mIHC_counts$subtypes_4, 
                                  col = list(HRD_status = c("high" = "salmon", "low/inter" = "seagreen1"),
                                             IM_status = c("0" = "#FFFFFF", "1" = "#000000"),
                                             PAM50_subtype = c("Basal" = "#FF5733", "Non Basal" = "#1E90FF", "Unclassified" = "#F5F5DC"),
                                             Refined_subtypes = c("BL1" = "#FF0000", "BL2" = "#FFA500", "LAR" = "#FFFF00", "M" = "#32CD32")))


annotations_right <- rowAnnotation(CD8_clusters=mIHC_counts$CD8_clusters,
                                   CD20_clusters=mIHC_counts$CD20_clusters,
                                   CD4_clusters=mIHC_counts$CD4_clusters,
                                   Tumour_clusters=mIHC_counts$Tum_clusters,
                                   Mixed_clusters=mIHC_counts$Mix_clusters,
                                   Homogeneous_tumour=mIHC_counts$HTum_clusters,
                                   Homogeneous_immune=mIHC_counts$HImm_clusters,
                                   Homogeneous_mixed=mIHC_counts$HMix_clusters)


Heatmap(log(mIHC_counts[,c("CD4", "CD8", "CD20", "CD4_FOXP3", "CD8_FOXP3", "PAN.CK", "Other")]+1),
        name= "Log cell counts",
        left_annotation = annotations_left,
        right_annotation = annotations_right,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        clustering_distance_rows = "euclidea")


# Combine data from the 2 different cores of the same sample

#
# COMPARING CORE TO CORE COUNTS AND GENERATING DATA POER SAMPLE
#

#Defining phenotypes
my_phenotypes <- c("CD8", "CD4", "CD20", "CD68", "CD4_FOXP3", "CD8_FOXP3", "PAN.CK", "Other")

#Getting sample ids
mIHC_counts$PDID <- annotation_file[mIHC_counts$TMArQ_CORE_ID, "PDid"]
my_samples <- unique(mIHC_counts$PDID)[!is.na(unique(mIHC_counts$PDID))]

# Generating dataframe to store the mean counts per sample
sample_counts_df <- data.frame(matrix(nrow=length(my_samples), ncol=length(my_phenotypes),  dimnames = list(my_samples, my_phenotypes)))

# Generating dataframe to store the correlation values between cores
correlation_dataframe <- as.data.frame(matrix(ncol = 2))
colnames(correlation_dataframe) <- c("Phenotype", "Correlation")


# Iterating through phenotypes and sample ids to fill in the matrices
for (phenotype in my_phenotypes) {
  
  df_of_counts <- as.data.frame(matrix(ncol = 2))
  
  for (sample in my_samples) {
    
    # Storing the counts of the cores belomnging to the same sample
    df_of_counts <- na.omit(rbind(df_of_counts, c(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1, phenotype] +1),
                      log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,phenotype])+1)))
    
    
    # Determining the mean counts for each phenotype and adding it to the dataframe
    #Using na.omit for the cases in whoch there is only a single core per sample
    sample_counts_df[sample, phenotype] <- mean(na.omit(c(mIHC_counts[which(mIHC_counts$PDID==sample),][1,phenotype]),
                                                          mIHC_counts[which(mIHC_counts$PDID==sample),][2,phenotype]))
    
  }
  
  correlation_dataframe <- rbind (correlation_dataframe, c(phenotype, cor(method="spearman", x=df_of_counts$V1, y=df_of_counts$V2)))
  
}

correlation_dataframe <- na.omit(correlation_dataframe)


#Plotting Spearman correlations
ggplot(data=correlation_dataframe, aes(x=as.factor(Phenotype), y=as.numeric(Correlation))) +
         geom_col() +
         ylim(0,1) +
         theme_classic() +
         labs(x="Cell phenotype", y="Core to core correlation (Spearman)") +
         theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
               text=element_text(size=12))


#Plotting scatterplots of each phenotype
plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")


for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD8"]+1),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD8"]+1),
         pch=16)


  
}

cd8_df <- na.omit(cd8_df)

plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")
for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD4"] +1 ),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD4"] +1),
         pch=16)
  
}



plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")
for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD20"] +1),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD20"] +1),
         pch=16)
  
}

plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")
for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD68"] +1),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD68"] +1),
         pch=16)
  
}

plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")
for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD4_FOXP3"] +1),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD4_FOXP3"] +1),
         pch=16)
  
}

plot(NA,xlim=c(0,10), ylim=c(0,10), xlab="Core 1 counts", ylab="Core 2 counts")
for (sample in my_samples) {
  
  points(log(mIHC_counts[which(mIHC_counts$PDID==sample),][1,"CD8_FOXP3"] +1),
         log(mIHC_counts[which(mIHC_counts$PDID==sample),][2,"CD8_FOXP3"] +1),
         pch=16)
  
}


#
# GENERATING HEATMAP FOR SAMPLES AND ANNOTATION (COMBINING CORES)
#

#Adding annotation to sample_counts dataframe
sample_counts_df$TILs <- as.numeric(sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"TILs"][1]))
sample_counts_df$ASCAT <- as.numeric(sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"ASCAT_TUM_FRAC"][1]))
sample_counts_df$HRD_st <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"HRD.2.status"][1])
sample_counts_df$PAM50 <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"PAM50_NCN"][1])
sample_counts_df$IM_status <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"TNBCtype_IMpositive"][1])
sample_counts_df$Lehman_4 <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"TNBCtype4_n235_notPreCentered"][1])
sample_counts_df$Lehman_7 <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"TNBCtype_org"][1])


# Defining annnotation columns for the heatmap
annotations <- rowAnnotation(TILs=sample_counts_df$TILs,
                             HRD_status=sample_counts_df$HRD_st,
                             PAM50_subtype=sample_counts_df$PAM50,
                             IM_status=sample_counts_df$IMstatus,
                             Refined_subtypes=sample_counts_df$Lehman_4)

Heatmap(log(sample_counts_df[,c("CD4", "CD8", "CD20", "CD4_FOXP3", "CD8_FOXP3", "PAN.CK", "Other")]+1),
        name= "Log cell counts",
        left_annotation = annotations,
        cluster_columns = FALSE,
        show_row_names = FALSE)


#
# ANALYSING SURVIVAL FROM SAMPLES (COMBINING CORE DATA)
#

#Adding annotation to sample_counts dataframe
sample_counts_df$time <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"time"][1])
sample_counts_df$status <- sapply(rownames(sample_counts_df), function(id) annotation_file[which(annotation_file$PDid==id),"event"][1])


# GENERATING FOREST PLOT

# Combine into a data frame
surv_data <- data.frame(time = as.numeric(sample_counts_df$time), 
                        status = sample_counts_df$status,
                        PAN.CK = sample_counts_df$PAN.CK,
                        CD8 = sample_counts_df$CD8,
                        CD4 = sample_counts_df$CD4,
                        CD20 = sample_counts_df$CD20,
                        CD68 = sample_counts_df$CD68,
                        CD4_FOXP3 = sample_counts_df$CD4_FOXP3,
                        CD8_FOXP3 = sample_counts_df$CD8_FOXP3,
                        Other = sample_counts_df$Other)

# Perform survival analysis
result = survivalAnalysis::analyse_multivariate(data = surv_data,
                                                time_status = c("time", "status"),
                                                covariates = c("CD8", "CD4", "CD20", "CD68", "CD4_FOXP3", "CD8_FOXP3"))

survivalAnalysis::forest_plot(result)


# GENERATING KAPLAND MEYER PLOTS

# Effect of T cyt
cutoff <- mean(sample_counts_df$CD8)
surv_data$Tcyt_status <- as.factor(sapply(surv_data$CD8, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Tcyt_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval = TRUE, ylab = "IDFS probability")

# Effect of T h
cutoff <- mean(sample_counts_df$CD4)
surv_data$Th_status <- as.factor(sapply(surv_data$CD4, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Th_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval = TRUE, ylab = "IDFS probability")

# Effect of B
cutoff <- mean(sample_counts_df$CD20)
surv_data$B_status <- as.factor(sapply(surv_data$CD20, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ B_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval = TRUE, ylab = "IDFS probability")

# Effect of Treg_cd4
cutoff <- mean(sample_counts_df$CD4_FOXP3)
surv_data$Treg_cd4_status <- as.factor(sapply(surv_data$CD4_FOXP3, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Treg_cd4_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval = TRUE, ylab = "IDFS probability")

# Effect of Treg_cd8
cutoff <- mean(sample_counts_df$CD8_FOXP3)
surv_data$Treg_cd8_status <- as.factor(sapply(surv_data$CD8_FOXP3, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Treg_cd8_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval=TRUE, ylab = "IDFS probability")

# Effect of CD68
cutoff <- mean(sample_counts_df$CD68)
surv_data$Mac_status <- as.factor(sapply(surv_data$CD68, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Mac_status, data = surv_data), risk.table = TRUE, conf.int = TRUE, pval = TRUE, ylab = "IDFS probability")


# Effect of all the immune cells
cutoff <- mean(sample_counts_df$CD8 + sample_counts_df$CD8_FOXP3 + sample_counts_df$CD4_FOXP3 + sample_counts_df$CD4 + sample_counts_df$CD20 + sample_counts_df$CD68)
surv_data$Immune_cell_count <- as.factor(sapply(surv_data$CD8_FOXP3 + surv_data$CD8 + surv_data$CD4_FOXP3 + surv_data$CD4 + surv_data$CD20 + surv_data$CD68, function(val) {if (val > cutoff) {"high"} else {"low"}}))
ggsurvplot(survival::survfit(Surv(time = time, event = status) ~ Immune_cell_count, data = surv_data), risk.table = TRUE, conf.int = TRUE, ylab = "IDFS probability", P_val=TRUE)



#
# PLOTTING SUBTYPES VS COUNTS
#


# Transforming dataframe to be plotted
sample_counts_df$sample_id <- rownames(sample_counts_df)

# Reassigning IM status to negative and positive
sample_counts_df$IM_status <- sapply(sample_counts_df$IM_status, function(val) {
  if (is.na(val)) {"NA"} 
  else if (val==1) {"Positive"}
  else {"Negative"}
})

# Gather the dataframe
sample_dataframe_long <- sample_counts_df %>% 
  gather(key = "Phenotype", value = "Count", c(1:8,))

# Defining order for plotting
lehman_7_order <- c("IM", "MSL", "BL1", "BL2", "LAR", "M", "UNS")
lehman_4_order <- c("BL1", "BL2", "LAR", "M")
im_status_order <- c("Negative", "Positive")
phenotype_order <- c("PAN.CK", "Other", "CD4", "CD8", "CD68", "CD20", "CD4_FOXP3", "CD8_FOXP3")

sample_dataframe_long$Lehman_7 <- factor(sample_dataframe_long$Lehman_7, levels = lehman_7_order)
sample_dataframe_long$Lehman_4 <- factor(sample_dataframe_long$Lehman_4, levels = lehman_4_order)
sample_dataframe_long$IM_status <- factor(sample_dataframe_long$IM_status, levels = im_status_order)
sample_dataframe_long$Phenotype <- factor(sample_dataframe_long$Phenotype, levels = phenotype_order)


# Lehman 7
ggplot(data=na.omit(sample_dataframe_long), aes(x=as.factor(Lehman_7), y=log(Count + 1), fill=as.factor(Phenotype))) +
  geom_boxplot() +
  labs(x="Lehman original subtypes", y="Log cell counts", fill="Cell phenotype") +
  theme_classic()

# Lehman 4
ggplot(data=na.omit(sample_dataframe_long), aes(x=as.factor(Lehman_4), y=Count, fill=as.factor(Phenotype))) +
  geom_boxplot() +
  labs(x="Lehman refined subtypes", y="Log cell counts", fill="Cell phenotype") +
  theme_classic()

# IM status
ggplot(data=na.omit(sample_dataframe_long), aes(x=as.factor(IM_status), y=log(Count + 1) , fill=as.factor(Phenotype))) +
  geom_boxplot() +
  labs(x="Lehman refined subtypes", y="Log cell counts", fill="Cell phenotype") +
  theme_classic()

# IM status
ggplot(data=na.omit(sample_counts_df), aes(x=as.factor(IM_status), y=log(CD8_FOXP3))) +
  geom_violin() +
  geom_boxplot(aes(width=0.6)) +
  labs(x="Lehman refined subtypes", y="Log cell counts", fill="Cell phenotype") +
  theme_classic()


#
# Plotting cibersort fractions VS groups
#


# Getting data of interest from annotation file
data_of_interest <- na.omit(annotation_file[,c("TNBCtype_IMpositive", "TNBCtype4_n235_notPreCentered", "CibersortX.Tcell", "CibersortX.Bcell", "CibersortX.macrophage")])
data_of_interest$TNBCtype_IMpositive <- sapply(data_of_interest$TNBCtype_IMpositive, function(val) {if (val==0) {"Negative"} else if (val==1) {"Positive"}})

data_of_interest <- gather(data_of_interest, key = "Cibersort", value = "Fraction", c(3,4,5))


ggplot(data_of_interest, aes(x=TNBCtype4_n235_notPreCentered, y=as.numeric(Fraction), fill=Cibersort)) +
  geom_boxplot() +
  labs(x="Lehman refined subtypes", y="Cibersort fraction") +
  theme_classic()


ggplot(data_of_interest, aes(x=as.factor(TNBCtype_IMpositive), y=as.numeric(Fraction), fill=Cibersort)) +
  geom_boxplot() +
  labs(x="IM status", y="Cibersort fraction") +
  theme_classic()

