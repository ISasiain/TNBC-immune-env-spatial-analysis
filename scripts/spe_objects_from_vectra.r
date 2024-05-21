#!/usr/bin/env Rscript

 # Description. This script generates spe objects from the phenotyped Vectra images.
 #
 # Usage. spe_objects_from_vectra.r -m [comma separated list] -i [path] -o [path]
 #
 # Arguments: -i [path] Path of the directory containing the coordinate files of each marker and sample.
 #            -m [comma separated list] Comma separated list of the markers to analyse.
 #            -o [path] Path of the directory to the file in which the output files will be stored.
 # Output: spe objects in the output directory.

#
# Loading required libraries
#

options(repos = "https://cran.r-project.org/")

#Installing packages if required
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")


if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")

suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(SPIAT))
suppressPackageStartupMessages(library(optparse))

#
# CONGIFIGURING COMMAND LINE OPTIONS
#

argument_list <- list(
  
  make_option(c("-m", "--markers"), type="character",  
              help="Comma separated list of the markers to analyse",
              metavar = "[comma separated list]"),
  
  make_option(c("-i", "--input_path"), type="character",
              help="Path of the directory containing the coordinate files of each marker and sample.",
              metavar = "[path]"),
  
  make_option(c("-o", "--output_path"), type="character",
              help="Path of the directory to the file in which the output files will be stored.",
              metavar = "[path]")
  
)


arguments <- parse_args(OptionParser(option_list=argument_list, 
                                     description="This program generates spe objects from the phenotyped Vectra images."))


#
# Generating file
#

files_vector <- list.files(arguments$input_path, pattern = "*intensity.csv")
my_markers <- strsplit(arguments$markers, ",")[[1]]

for (file_path in files_vector) {
  
  my_example_file <- read.csv(paste0(arguments$input_path, file_path))

  #Splitted filename
  splitted_filename <- strsplit(file_path, "_")[[1]]
  
  # Getting information to fill the columns of teh dataframe
  
  # Getting the core's id
  coreinfo <- strsplit((strsplit(splitted_filename[3], '[][]+')[[1]][2]), ",")[[1]]
  
  my_block <- strsplit(splitted_filename[2], "")[[1]][1]
  my_core <- paste0(coreinfo[3], coreinfo[2])


  core_name <- paste0("Core_", my_block, "_", my_core)
  
  cell_ids <- sapply(rownames(my_example_file), function(id) paste0("Cell_", id))
  intensity_matrix <- t(as.matrix(my_example_file[,my_markers]))
  colnames(intensity_matrix) <- cell_ids
  
  
  #Generate spe object
  my_spe <- format_image_to_spe(
    format="general",
    intensity_matrix = intensity_matrix,
    phenotypes = my_example_file$Phenotype,
    coord_x = my_example_file$X_coor,
    coord_y = my_example_file$Y_coor
  )
  
  # Saving file
  
  filename <- paste0(arguments$output_path, core_name, "_spe.object.rds")
  
  saveRDS(my_spe, file = filename)
  
}

cat ("\n\nAll files have been generated successfully.\n")
