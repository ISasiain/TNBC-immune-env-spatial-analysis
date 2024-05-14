#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexHeatmap)
library(SPIAT)
library(SpatialExperiment)
library(ComplexHeatmap)


#
# GETTING PATHS TO SPE OBJECTS
#

# Defining path to spe_objects
path_to_spes <- "/media/isc/Data1/spe_objects/"

#Getting paths to spe files
spe_files <- list.files(path_to_spes)

#
# READING ANNOTATION FILE
#

annotation_file <- read.csv2("/home/isc/Spatial_immune_env/data_from_suze/data/supplData_withimages.csv")
rownames(annotation_file) <- annotation_file$uid

#
# READING CELL COUNTS 
#

mIHC_counts1 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_1_dataframe.csv")
mIHC_counts2 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_2_dataframe.csv")
mIHC_counts3 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_3_dataframe.csv")
mIHC_counts4 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_4_dataframe.csv")
mIHC_counts5 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_5_dataframe.csv")


mIHC_counts <- rbind(mIHC_counts1, mIHC_counts2, mIHC_counts3, mIHC_counts4, mIHC_counts5)
rownames(mIHC_counts) <- mIHC_counts$CORE_ID

#
# DETERMININING PARAMETERS: Average minimum distance
#


# Defining phenotype list
phenotypes <- c("CD8", "CD4", "CD20", "CD68", "CD4_FOXP3", "CD8_FOXP3", "PAN-CK", "Other")

# Generating dataframe to store distances
distances_df <- as.data.frame(matrix(spe_files, nrow=length(spe_files), ncol=1))
colnames(distances_df) <- c("Core")

# Generating average minimum distance between phenotypes in each core
for (phenotype1 in phenotypes) {
  for (phenotype2 in phenotypes[!phenotype1 ==  phenotypes]) {
    
    #Generating colnames of the data frame
    colname <- paste0(phenotype1, "/", phenotype2)
    
    #Appending column to data frame
    distances_df[[colname]] <- 
      sapply(spe_files, function(file) {
        print(file)
        
        tryCatch({
          result <- calculate_summary_distances_between_celltypes(
            calculate_minimum_distances_between_celltypes(
              spe_object = readRDS(paste0(path_to_spes, file)),
              cell_types_of_interest = c(phenotype1, phenotype2),
              feature_colname = "Phenotype"
            )
          )[1, 5]
          
          return(result)
        }, error = function(e) {
          return(NaN)
        })
      })
    
  }
}


#
# GETTING BLOCK ID AND TMARQ BLOCK ID
#

# Generating Block ID 
distances_df$Core_ID <- sapply(distances_df$Core, 
       function(name) {
         paste0(
           "BLOCK_",
           strsplit(name, "_")[[1]][2],
           "|",
           strsplit(name, "_")[[1]][3]
         )
      })

# Adding TMArQ ID
distances_df$TMArQ_ID <- sapply(distances_df$Core_ID, 
                               function(id) {
                                  mIHC_counts[id, "TMArQ_CORE_ID"]
                               })


#
# ADDING ANNOTATIONS TO DATA FRAME
#

distances_df$Lehman_7 <- annotation_file[distances_df$TMArQ_ID, "TNBCtype_org"]
distances_df$IM_status <- annotation_file[distances_df$TMArQ_ID, "TNBCtype_IMpositive"]
distances_df$Lehman_4 <- annotation_file[distances_df$TMArQ_ID, "TNBCtype4_n235_notPreCentered"]


#
# SUMMARISING DISTANCES PER GROUP. MEDIAN
#

# Determining phenotype combinations
phenotype_combinations <- c()

for (phenotype1 in phenotypes) {
  for (phenotype2 in phenotypes[!phenotype1 ==  phenotypes]) {
    
    phenotype_combinations <- c(phenotype_combinations,
                                paste0(phenotype1, "/", phenotype2))
    
  }
}

# Generate dataframe to store the data
distances_summary <- as.data.frame(matrix(nrow=7, ncol=length(phenotype_combinations)))
colnames(distances_summary) <- phenotype_combinations
rownames(distances_summary) <- c("All samples", 
                                 "IM +", "IM -",
                                 "BL1", "BL2", "LAR", "Refined Lehman: M")

# Original_values
distances_summary["All samples",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[,combination], na.rm = TRUE)
  
})


# IM positive 
distances_summary["IM +",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$IM_status == 1),combination], na.rm = TRUE)
  
})

# IM negative
distances_summary["IM -",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$IM_status == 0),combination], na.rm = TRUE)
  
})

# L4_BL1
distances_summary["BL1",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$Lehman_4 == "BL1"),combination], na.rm = TRUE)
  
})

# L4_BL2
distances_summary["BL2",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$Lehman_4 == "BL2"),combination], na.rm = TRUE)
  
})

# L4_LAR
distances_summary["LAR",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$Lehman_4 == "LAR"),combination], na.rm = TRUE)
  
})

# L4_M
distances_summary["M",] <- sapply(phenotype_combinations, function(combination) {
  
  median(distances_df[which(distances_df$Lehman_4 == "M"),combination], na.rm = TRUE)
  
})

# # L7_IM
# distances_summary["Original Lehman: IM",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "IM"),combination], na.rm = TRUE)
#   
# })
# 
# # L7_MSL
# distances_summary["Original Lehman: MSL",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "MSL"),combination], na.rm = TRUE)
#   
# })
# 
# # L7_BL1
# distances_summary["Original Lehman: BL1",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "BL1"),combination], na.rm = TRUE)
#   
# })
# 
# # L7_BL2
# distances_summary["Original Lehman: BL2",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "BL2"),combination], na.rm = TRUE)
#   
# })
# 
# # L7_LAR
# distances_summary["Original Lehman: LAR",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "LAR"),combination], na.rm = TRUE)
#   
# })
# 
# 
# # L7_M
# distances_summary["Original Lehman: M",] <- sapply(phenotype_combinations, function(combination) {
# 
#   median(distances_df[which(distances_df$Lehman_7 == "M"),combination], na.rm = TRUE)
#   
# })
# 
# 
# # L7_UNS
# distances_summary["Original Lehman: UNS",] <- sapply(phenotype_combinations, function(combination) {
#   
#   median(distances_df[which(distances_df$Lehman_7 == "UNS"),combination], na.rm = TRUE)
#   
# })


#
# PLOTTING THE OBTAINED RESULTS
#

# Defining row groups
grouping_rows <- c("IM status", "IM status", "Refined subtypes", "Refined subtypes", "Refined subtypes", "Refined subtypes")
grouping_cols <- factor(c(rep("CD8",7), rep("CD4",7), rep("CD20",7), rep("CD68",7), rep("CD4_FOXP3",7), rep("CD8_FOXP3",7), rep("PAN-CK",7), rep("Other",7)), levels = c("CD8", "CD4", "CD20", "CD68", "CD4_FOXP3", "CD8_FOXP3", "PAN-CK", "Other"))

Heatmap(scale(distances_summary[rownames(distances_summary)[2:length(rownames(distances_summary))],], center = distances_summary["All samples",]),
        name ="Scaled median AMD",
        row_order = rownames(distances_summary)[2:length(rownames(distances_summary))],
        column_order = colnames(distances_summary),
        row_split = grouping_rows,
        column_split = grouping_cols)
