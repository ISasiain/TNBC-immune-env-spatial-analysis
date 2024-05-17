#!/usr/bin/Rscript

#
# DOWNLOADING AND LOADING REQUIRED PACKAGES
#

options(repos = "https://cran.r-project.org/")

if (!require("survival", quietly = TRUE))
  install.packages("survivalr")
library(survival)

if (!require("survminer", quietly = TRUE))
  install.packages("survminer")
library(survminer)

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("ComplexHeatmap", quietly = TRUE))
  install.packages("ComplexHeatmap")
library(ComplexHeatmap)

if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
library(stringr)


#
# READING FILES
#

dins <- readRDS("/home/isc/Spatial_immune_env/Spatial_immune_env_analysis/detected_clusters/r100_DIN_sIHC.rds")
clustered_cells <- readRDS("/home/isc/Spatial_immune_env/Spatial_immune_env_analysis/detected_clusters/r100_DIN_sIHC_clusters_clustered_cells.rds")
core_clusters <- readRDS("/home/isc/Spatial_immune_env/Spatial_immune_env_analysis/detected_clusters/r100_DIN_sIHC_clusters.rds")
annotation_file <- readRDS("/home/isc/Spatial_immune_env/data_from_suze/data/supplData_withimages.csv")
rownames(annotation) <- annotation$uid

# 
# GETTINBG BLOCK IDS FROM PATHS TO USE AS NAMES
#

names(clustered_cells) <- str_extract(names(clustered_cells), "BLOCK_\\d+\\|\\w\\d+")
names(core_clusters) <- str_extract(names(core_clusters), "BLOCK_\\d+\\|\\w\\d+")
names(dins) <- str_extract(names(dins), "BLOCK_\\d+\\|\\w\\d+")


#
# PARSING FILES. GETTING DISMISSED, HOMOGENEOUS AND HETEROGENEOUS CORES
#

# Getting cores with low density
dismissed_cores = c()

# Iterating over the samples
for (element in names(core_clusters)) {
  #Checking if the length of the list's element is 1
  if (length(core_clusters[[element]]) <= 1) {
    #Appending element to vector
    dismissed_cores <- c(dismissed_cores, element)
  }
}


# Getting homogeneous cores
homogeneous_cores = c()

# Iterating over the samples
for (element in names(core_clusters)) {
  #Checking if the length of the list's element is 1
  if (length(core_clusters[[element]]) == 2) {
    if (core_clusters[[element]]$"type" == "Homogeneous") {
      #Appending element to vector
      homogeneous_cores <- c(homogeneous_cores, element)
    }
  }
}


# Getting non-homogeneous cores. All the core IDs not in the previous vectors
non_homogeneous_cores = names(core_clusters)[!names(core_clusters) %in% c(homogeneous_cores, dismissed_cores)]  


#
# ANALYSING COMPOSITION OF NON-HOMOGENEOUS CORES
#

# MDefine markers generate epty dataframe to store the data
markers <- c("p53", "CD3", "CD8", "CD20", "CD4", "CKPAN", "FOXP3", "H2AXp", "CD68")
cluster_dens_df <- data.frame()


# Getting only the din matrices of the non-homogeneous cores

for (core in non_homogeneous_cores) {
  
  din_matrix <- dins[[core]]
  clusters <- clustered_cells[[core]]
  cluster_annotation <- core_clusters[[core]]
  
  
  if (length(cluster_annotation[["clusters"]]) != 0) {
    for (cluster in 1:length(cluster_annotation)) {
    
      new_row <- c(
        core,
        cluster,
        cluster_annotation[["clusters"]][[cluster]][["type"]],
        sapply(markers, function(marker) {
          
          #Using tryCatch to prevent the script from failing if one of the markers is not included 
          #in one of the DIN matrices
          tryCatch({
            value <- median(din_matrix[clusters$"cluster" == cluster, marker])
            ifelse(is.null(value), 0, value)
          }, error = function(e) {
            0
          })
          
        })
      )
    
    
      print(new_row)
      
      cluster_dens_df <- rbind(new_row, cluster_dens_df)
      
    }
    
  }
  
}


# Adding colnames to the cluster_dens_df
colnames(cluster_dens_df) <- c("Sample", "Cluster", "Type", markers)


# Generating heatmap of scaled median density


#Tranforming values to numeric and factors
cluster_dens_df$p53 <- as.numeric(cluster_dens_df$p53)
cluster_dens_df$CD3 <- as.numeric(cluster_dens_df$CD3)
cluster_dens_df$CD8 <- as.numeric(cluster_dens_df$CD8)
cluster_dens_df$CD20 <- as.numeric(cluster_dens_df$CD20)
cluster_dens_df$CD68 <- as.numeric(cluster_dens_df$CD68)
cluster_dens_df$CD4 <- as.numeric(cluster_dens_df$CD4)
cluster_dens_df$CKPAN <- as.numeric(cluster_dens_df$CKPAN)
cluster_dens_df$FOXP3 <- as.numeric(cluster_dens_df$FOXP3)
cluster_dens_df$H2AXp <- as.numeric(cluster_dens_df$H2AXp)

cluster_dens_df$Type <- as.factor(cluster_dens_df$Type)


# Generating scatterplots of the markers used for clustering vs p53

ggplot(data=cluster_dens_df, aes(x=p53, y=CD3)) +
  geom_point(aes(color=Type)) +
  xlab(paste0("Median DIN of p53")) +
  ylab(paste0("Median DIN of CD3")) +
  theme_classic()

ggplot(data=cluster_dens_df, aes(x=p53, y=CD4)) +
  geom_point(aes(color=Type)) +
  xlab(paste0("Median DIN of p53")) +
  ylab(paste0("Median DIN of CD4")) +
  theme_classic()

ggplot(data=cluster_dens_df, aes(x=p53, y=CD8)) +
  geom_point(aes(color=Type)) +
  xlab(paste0("Median DIN of p53")) +
  ylab(paste0("Median DIN of CD8")) +
  theme_classic()

ggplot(data=cluster_dens_df, aes(x=p53, y=CD20)) +
  geom_point(aes(color=Type)) +
  xlab(paste0("Median DIN of p53")) +
  ylab(paste0("Median DIN of CD20")) +
  theme_classic()


#Generating sorting vector
sorting_vec <- order(cluster_dens_df[,"Type"])

#Scaling median DINs
cluster_dens_df[,seq(4,12)] <- scale(cluster_dens_df[,seq(4,12)])

my_annotation <- rowAnnotation(
  Cluster_type = cluster_dens_df[,"Type"],
  col = list(Cluster_type = c("TUMOUR" = "magenta4", "IMMUNE" = "lightsalmon", "MIXED" = "beige"))
)

Heatmap(
  cluster_dens_df[, markers],
  cluster_rows = FALSE,
  row_order = sorting_vec,
  row_title = "Total clusters detected",
  column_title = "Markers",
  name = "Scaled median DIN",
  show_row_names = FALSE,
  column_order = markers,
  left_annotation = my_annotation
)

