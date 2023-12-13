rm(list = ls())
setwd("/home/isc/Spatial_immune_env")
#Installing packages if requitred
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("stringr", quietly = TRUE))
  install.packages("stringr")

if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")

if (!require("R.utils", quietly = TRUE))
  install.packages("R.utils")

if (!require("alphahull", quietly = TRUE))
  install.packages("alphahull")

if (!require("plotly", quietly = TRUE))
  install.packages("plotly")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("lme4", quietly = TRUE))
  install.packages("lme4")


#Loading packages if required
library(SpatialExperiment)
library(SPIAT)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(alphahull)
library(plotly)
library(tidyr)
library(lme4)


#
# DEFINING PARAMETERS 
#

#Generate a list to append summary statistics of each image
list_of_summary_statistics <- list()

#Generate a list of the paths to the image objects
path_to_files <- "data_from_suze/data/spiat/"
spe_objects <- list.files(path_to_files)

#Defining reference and target markers

markers1 <- c("p53", "CD3", "CD8", "CD20")
markers2 <- c("p53", "CD3", "CD8", "CD20")

for (marker_ref in markers1) {
  
  for (marker_target in markers2[markers2 != marker_ref]) {
    
    print(marker_ref)
    print(marker_target)
       
    reference_markers <- c(marker_ref)
    target_markers <- c(marker_target)
    
    
    #
    # GENERATING METRICS
    #
    
      for (spe in spe_objects) {
        
        sample <- readRDS(paste(path_to_files, spe, sep=""))
        
        #Obtaining ID to use as identifier on the list
        sample_ID <- paste("BLOCK_", strsplit(spe, "_")[[1]][3], sep="")
        
        #Only calculate parameters when the numbers of cells with a certain phenotype
        #is higher than 2.
        if (reference_markers %in% sample$Phenotype & 
            target_markers %in% sample$Phenotype &
            table(sample$Phenotype)[reference_markers] > 2 &
            table(sample$Phenotype)[target_markers] > 2) {
          
          
          #Estimate summary statistics of each sample
          list_of_summary_statistics[[sample_ID]] <- list(
            
            #Determining phenotype counts of the sample
            "Cell_Counts" = 
              table(sample$Phenotype), 
            
            #Determining phenotype proportions on the sample
            "Cell_Proportions" = 
              100*table(sample$Phenotype)/length(sample$Phenotype),
            
            # Average minimum distance between reference and targets
            "AMD" = 
              calculate_summary_distances_between_celltypes(
                calculate_minimum_distances_between_celltypes(
                  sample,
                  cell_types_of_interest = c(reference_markers, target_markers), 
                  feature_colname = "Cell.Type")
              ),
            
            # Cells in the neighborhood. Only with a single marker
            "CIN" = average_percentage_of_cells_within_radius(
              sample, 
              reference_celltype = reference_markers, 
              target_celltype = target_markers, 
              feature_colname = "Cell.Type", 
              radius = 100),
            
            # Area under curve of the K cross function
            "AUC" = calculate_cross_functions(sample, 
                                              method = "Kcross", 
                                              cell_types_of_interest = c(reference_markers, target_markers), 
                                              feature_colname = "Cell.Type"),
            
            # Entropy (attraction measure)
            "Entropy" = entropy_gradient_aggregated(
                          spe_object = sample,
                          radii = c(50, 100, 150, 200, 250, 300), 
                          cell_types_of_interest = c(reference_markers, target_markers),
                          feature_colname = "Cell.Type")
            
          )
          
        } else {
          
          #Estimate summary statistics of each sample
          list_of_summary_statistics[[sample_ID]] <- list(NA)
          
        }
        
        saveRDS(list_of_summary_statistics, paste0(reference_markers,"_", target_markers,"_spatial_stats.RData"))
        
      }
    }
  }
  

