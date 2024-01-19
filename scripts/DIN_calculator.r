#!/usr/bin/Rscript

#
# LOADING THE REQUIRED PACKAGES
#

# Set a specific CRAN mirror
options(repos = "https://cran.r-project.org/")

if (!require("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse)

if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")
library(SpatialExperiment)

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")
library(SPIAT)

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("stats", quietly = TRUE))
  install.packages("stats")
library(stats)

if (!require("fpc", quietly = TRUE))
  install.packages("fpc")
library(fpc)

if (!require("foreach", quietly = TRUE))
  install.packages("foreach")
library(foreach)

if (!require("doParallel", quietly = TRUE))
  install.packages("doParallel")
library(doParallel)


#
# CONFIGURE COMMAND LINE ARGUMENTS
#

argument_list <- list(

    make_option(c("-o", "--spe_objects"), 
                type="character", 
                help="Comma separated spe objects' directories.",
                metavar="[COMMA_SEPARATED_PATHS]"),
    
    make_option(c("-c", "--cores"), 
              type="integer", 
              default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[NUMBER]"))

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program generates DIN values of the analysed samples."))

#
# CONFIGURATE PARALLELIZATION
#

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores)  
registerDoParallel(cl)  


#
# READING SPE OBJECTS
#

#Splitting command line argument to generate a vector of the spe objects's paths

print("Reading SPE objects...")

spe_names <- strsplit(arguments$spe_objects, split=",")[[1]]


#Storing spe objects in a list. Using parLapply to parallelize the process

file_name_list <- foreach(file=spe_names) %dopar% {

    list(file, 
         readRDS(file))

}

#Rearrange list
file_list <- lapply(file_name_list, function(item) item[[2]])
names(file_list) <- sapply(file_name_list, function(item) item[[1]])

#
# GENERATING CELL DENSITY VALUES
#


#Generating function to calculate density values
calculate_density_matrix <- function(spe_object, radius) {

            #Getting all the markers included in the spe object
            markers <- unique(spe_object$"Cell.Type")


            #Generating matrix to store the denisty in the neighborhood (DIN) values
            number_of_columns <- 2 + length(unique(markers))
            number_of_rows <- length(spe_object$"Cell.ID")


            DIN_matrix <- matrix(NA, nrow=number_of_rows, ncol=number_of_columns)

            rownames(DIN_matrix) <- spe_object$"Cell.ID"
            colnames(DIN_matrix) <- c("X_coor", "Y_coor", unique(spe_object$"Cell.Type"))

            #Generating data frame with all the cellIDs, Cell.Type and coordinates
            spatial_data <- cbind(spe_object@colData[, c("Cell.ID", "Cell.Type")], spatialCoords(spe_object))


            #Getting spatial coordinates
            DIN_matrix[,c("X_coor",  "Y_coor")] <- spatialCoords(spe_object)
            
            #Removing cells without defined coordinates. It happens in rare cases
            spatial_data <- na.omit(spatial_data)
            DIN_matrix <- DIN_matrix[spatial_data$Cell.ID,]

            #Calculating DIN values
            for (cell_type in unique(markers)) {
            
            #Use sapply to generate a vector with the DINs of each cell type
            DIN_matrix[,cell_type] <- sapply(
                
                #Calculating DIN for each cell in the sample
                rownames(DIN_matrix), 
                
                function(cell) {
                
                #Defining parameters
                X_coor <- DIN_matrix[cell, "X_coor"]
                Y_coor <- DIN_matrix[cell, "Y_coor"]
                
                #Creating or reassigning a column in spatial_data dataframe to calculate
                #the euclidean distance between cells
                spatial_data$"Distance" <-
                    (spatial_data$Cell.X.Position - X_coor)^2 + 
                    (spatial_data$Cell.Y.Position - Y_coor)^2
              

                #Filtering the cells not included in the area of interest and that do not
                #belong to the category of interest. Counting the cells that match the conditions.
                return(nrow(spatial_data[spatial_data$"Distance" <= radius^2 & spatial_data$"Cell.Type" == cell_type,])/(pi*radius^2))
                
                }
            )

            }

    return(DIN_matrix)
}


print("Calculating density in the neighborhood matrices...")


density_matrix_list <- foreach(s=spe_names, .packages=c("SPIAT")) %dopar% {
  
  # Generating list with sample id and dnsity matrix
  list(strsplit(strsplit(s, "/")[[1]][4], "_for.spiat")[[1]][1],
       calculate_density_matrix(file_list[[s]], radius=100)
      )

}

# Stop the defined clusters
stopCluster(cl)

#Rearrange list
output_list <- lapply(density_matrix_list, function(item) item[[2]])
names(output_list) <- sapply(density_matrix_list, function(item) item[[1]])

#
# SAVING OUTPUT FILES
#

print("Saving output file...")

# Saving output file

saveRDS(output_list, file="all_samples_DIN.rds")

# Print message to show the end of the execution
cat("\n\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")