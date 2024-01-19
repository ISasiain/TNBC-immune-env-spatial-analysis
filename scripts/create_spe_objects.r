#!/usr/bin/Rscript

#
# INSTALL AND LOAD PACKAGES IF REQUIRED
#

options(repos = "https://cran.r-project.org/")

#Installing packages if required
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")


#Loading packages if required
suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(SPIAT))
suppressPackageStartupMessages(library(optparse))


#
# CONGIFIGURING COMMAND LINE OPTIONS
#

argument_list <- list(
  
  make_option(c("-c", "--cores"), type="integer", default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[number]"),
  
  make_option(c("-m", "--markers"), type="character",  
              help="Comma separated list of the markers to analyse",
              metavar = "[comma separated list]"),
  
  make_option(c("-a", "--annotation"), type="character", 
              help="Path to the annotation file paths containing the files containing each marker's coordinate. Include all the backslashes.",
              metavar = "[path]"),
  
  make_option(c("-p", "--path"), type="character",
              help="Path of the directory containing the coordinate files of each marker and sample.",
              metavar = "[path]")
  
)


arguments <- parse_args(OptionParser(option_list=argument_list, 
                                     description="This program generates spe object integrating files containing coordenates of different markers."))


#
# GENERATE SPE OBJECTS INTEGRATING DATA FROM DIFFERENT MARKERS
#

#Reading annotation file
annotation <- read.csv2(arguments$annotation, sep="\t")

#Getting vector of markers
markers <- strsplit(arguments$markers, ",")[[1]]

#Iterate through every sample, every row of the annotation file
for (row in 1:nrow(annotation)) {
  
  #Initializing data frame to store cell, coordinates and phenotype
  df_of_cells <- data.frame(matrix(ncol=3,nrow=0))
  
  for (marker in markers) {
    
    #Reading the file containing the marker coordinates from the annotation file and opening it
    col <- paste0(marker, "_image")
    coor_file <- read.csv2(paste0(arguments$path, annotation[row, col], "_coordinates.txt"), sep="\t")
    
    #Append all the data to df_of_cells dataframe
    df_of_cells <- rbind(df_of_cells, data.frame(coor_file$x, coor_file$y, rep(marker, nrow(coor_file))))
    
  }
  
  #Rename columns of the df_of_cells
  colnames(df_of_cells) <- c("X", "Y", "Phenotype")
  
  #Create cell ids and use them as rownames
  rownames(df_of_cells) <- sapply(1:nrow(df_of_cells), function(num) paste0("Cell_", num))



#
# GENERATING SPE OBJECT
#

# Generating to store the new SPE object as an rdata file.
pdid <- annotation[row, "PDid"]
uid <- annotation[row, "uid"]
name <- paste0(pdid, "_", uid, "_speGG.rds")

#Generating the intensity matrix for all the markers analysed. Using 20 as default value.
#Using the markers as rownames and the cell ids as colnames.
intensity_matrix <- matrix(NA, nrow=length(markers), ncol=nrow(df_of_cells))
rownames(intensity_matrix) <- markers
colnames(intensity_matrix) <- rownames(df_of_cells)

# Iterate through the markers analysed
for (marker in rownames(intensity_matrix)) {
  
  #Generate a vector containing the intensities of each marker in all the cells.
  intensity_matrix[marker,] <- sapply(colnames(intensity_matrix),
                                  FUN=function(cell) {
                                    #If the phenotype of the cell is the marker analysed
                                    if (marker==df_of_cells[cell, "Phenotype"]) {
                                      20 #Use the default intensity value
                                    } else {
                                      0 # Use 0 as intensity if it is not
                                    }
                                  }
                                )
                              }

#Creating SPE object using the generated data
my_spe <- format_image_to_spe(
          format = "general",
          intensity_matrix = intensity_matrix,
          phenotypes = df_of_cells$Phenotype,
          coord_x = df_of_cells$X,
          coord_y = df_of_cells$Y
        )


#Adding cell type to the generated spe object
my_spe <- define_celltypes(
  my_spe, 
  categories = markers,
  category_colname = "Phenotype", 
  names = markers,
  new_colname = "Cell.Type"
)

#Saving the generated spe object as an robject.
}