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

if (!require("diagram", quietly = TRUE))
  install.packages("diagram")


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
library(diagram)

#
# Defining functions to generat new metrics
#

cell_density_in_neighborhood <- function(spe_object, reference_cell, target_cell, radius) {
  
  # Combine cell ids, phenotypes, and coordinates of cells in a single data frame
  spatial_data <- cbind(J10_B1@colData[, c("Cell.ID", "Phenotype")], spatialCoords(J10_B1))
  
  # Filter the data frame to keep only the cells belonging to the reference and target phenotypes 
  reference <- spatial_data[spatial_data$Phenotype == reference_cell, ]
  target <- spatial_data[spatial_data$Phenotype == target_cell, ]
  
  # Raise an error if the phenotypes selected do not exist in the data   
  if (nrow(reference) == 0) {stop("The reference phenotype selected does not exist.")}
  if (nrow(target) == 0) {stop("The target phenotype selected does not exist.")}
  
  # Calculate number of cells in the neighborhood and append them to a named vector
  ncells_in_n <- sapply(1:nrow(reference), function(i) {
    
    ref_cell <- reference[i, ]  # Access each reference cell row-wise
    
    target$"Distance" <-
      (target$Cell.X.Position - ref_cell$Cell.X.Position)^2 + 
      (target$Cell.Y.Position - ref_cell$Cell.Y.Position)^2
    
    #Returning the points within the selected radius
    return(sum(target$"Distance" <= radius^2))
    
  })
  
  #Eliminate 1 cell (the reference) if the target and reference cells are the same
  if (target_cell == reference_cell) {ncells_in_n <- ncells_in_n -1}
  
  #Add names to vector
  names(ncells_in_n) <- reference$Cell.ID
  return(ncells_in_n / (pi*radius^2))
}

#
# Generating metrics
#

#Generate a list to append summary statistics of each image
list_of_summary_statistics <- list()

#Generate a list of the paths to the image objects
path_to_files <- "data_from_suze/data/spiat/"
spe_objects <- list.files(path_to_files)
#markers_to_analyse <- c("CD3", "CD8", "CD20", "p53")
markers_to_analyse <- c("p53")

for (marker_to_analyse in markers_to_analyse) {

  for (spe in spe_objects) {
    
    sample <- readRDS(paste(path_to_files, spe, sep=""))
    
    #Obtaining ID to use as identifier on the list
    sample_ID <- paste("BLOCK_", strsplit(spe, "_")[[1]][3], sep="")
    
    #Only calculate parameters when the numbers of cells with a certain phenotype
    #is higher than 2.
    if (marker_to_analyse %in% sample$Phenotype & table(sample$Phenotype)[marker_to_analyse] > 2) {
      
    #Estimate summary statistics of each sample
    list_of_summary_statistics[[sample_ID]] <- list(
      
      #Determining phenotype counts of the sample
      "Cell_Counts" = 
        table(sample$Phenotype), 
      
      #Determining phenotype proportions on the sample
      "Cell_Proportions" = 
        100*table(sample$Phenotype)/length(sample$Phenotype), 
      
      
      # Determining the cell of interest as reference. The radius was set to 200
      "CIN"= average_percentage_of_cells_within_radius(spe_object = sample, 
                                                       feature_colname = "Cell.Type", 
                                                       reference_celltype = marker_to_analyse, 
                                                       target_celltype = marker_to_analyse, 
                                                       radius = 100),
      
      # Average pairwise distance
      #"APD"= calculate_summary_distances_between_celltypes(
      #  calculate_pairwise_distances_between_celltypes(spe_object = sample, 
      #                                                 cell_types_of_interest = marker_to_analyse, 
      #                                                 feature_colname = "Cell.Type")
      #       ),
      
      # ANNI index
      "ANNI"= average_nearest_neighbor_index(spe_object = sample,
                                             reference_celltypes = marker_to_analyse, 
                                             feature_colname = "Cell.Type")
      
    )
    
    } else {
      
      #Estimate summary statistics of each sample
      list_of_summary_statistics[[sample_ID]] <- list(NA)
    
  }
  
  saveRDS(list_of_summary_statistics, paste0(marker_to_analyse, "_props_and_CIN.RData"))
  
  }
}
