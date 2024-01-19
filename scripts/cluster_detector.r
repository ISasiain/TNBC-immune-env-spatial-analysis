#!/usr/bin/Rscript

#
# DOWNLOADING AND LOADING REQUIRED PACKAGES
#

options(repos = "https://cran.r-project.org/")

if (!require("foreach", quietly = TRUE))
  install.packages("foreach")
library(foreach)

if (!require("doParallel", quietly = TRUE))
  install.packages("doParallel")
library(doParallel)

if (!require("NbClust", quietly = TRUE))
  install.packages("NbClust")
library(NbClust)

if (!require("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse)

if (!require("stats", quietly = TRUE))
  install.packages("stats")
library(stats)

#
# PARSING COMMAND LINE ARGUMENTS
#

argument_list <- list(

    make_option(c("-d", "--DIN_matrix_list"), 
                type="character", 
                help="The path to the .rds file containing the list of DIN matrices of the analysed 
                        samples must be entered here.",
                metavar="[PATH]"),

    make_option(c("-m", "--markers_to_analyse"), 
                type="character", 
                help="A comma separated list of the markers to be used to determine clusters should be entered here.",
                metavar="[PATH]"),
    
    make_option(c("-c", "--cores"), 
              type="integer", 
              default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[NUMBER]"))

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program generates cell clusters based on DIN values."))


#
# READING INPUT ARGUMENTS
#

# Loading input file
DIN_matrices <- readRDS(arguments$DIN_matrix_list)

# Generating vector of markers to analyse
markers <- strsplit(arguments$markers_to_analyse, ",")[[1]]

#
# CONFIGURATING PARALLELIZATION
#

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores)  
registerDoParallel(cl)  

#
# DETERMINING OPTIMAL NUMBER OF CLUSTERS
#

#Determining clusters

# Determining clusters

optimal_number_of_clusters <- foreach(sample=names(DIN_matrices), .packages="NbClust") %dopar% {
  
  #If the density of some of the markers to be used for clustering is not defined
  #the sample will be ignored.
  if (all(markers %in% colnames(DIN_matrices[[sample]]))) {
  
      # Generate optimal number of clusters based on the specified markers
    list(sample,
        length(unique(NbClust(
            DIN_matrices[[sample]][,markers],
            method="kmeans" #Using kmeans clustering
        )$Best.partition))
    )

  }

}

#Rearraging list
num_clusters_list <- lapply(optimal_number_of_clusters, function(item) item[[2]])
names(num_clusters_list) <- sapply(optimal_number_of_clusters, function(item) item[[1]])

print(num_clusters_list)

# Print warning message about the samples not used for clustering
not_clustered <- names(DIN_matrices)[!names(DIN_matrices) %in% names(num_clusters_list)]

if (length(not_clustered) != 0) {

  cat("\n\n**** WARNING *****\n")
  cat("\nThe following samples were not used for clustering as data for some of the markers was unavailable\n")
  print(not_clustered)

}


 output_clusters_list <- list()
 for (sample in names(num_clusters_list)) {

   output_clusters_list[[sample]] <- kmeans(DIN_matrices[[sample]][,markers], centers=num_clusters_list[[sample]])

 }

#Generate a cluster list base on the optimal number of clusters detected
#clusters_list <- foreach(sample=names(num_clusters_list), .packages="stats") %do% {
#
#    list(sample,
#      kmeans(DIN_matrices[[sample]][,markers], centers=num_clusters_list[[sample]]),
#      )
#}

# Rearranging list
#output_clusters_list <- lapply(clusters_list, function(item) item[[2]])
#names(output_clusters_list) <- sapply(clusters_list, function(item) item[[1]])


#
# CLUSTER EVALUATION
#

# # Generate a function to evaluate the detected clusters
# cluster_evaluation <- function(DIN_matrix, clusters_object, markers) {

#   # Obtain the number of clusters
#   cluster_num <- unique(unname(clusters_object))

#   # Generate data frame to store the comparison between clusters detected
#   cluster_comp <- data.frame()
#   colnames(cluster_comp) <- cluster_num
#   rownames(cluster_comp) <- cluster_num

#   # Compare clusters
#   for(my_row in cluster_num) {
    
#     #Filling th rows of the cluster_comp data frame
#     cluster_comp[cluster,] <- sapply(
#       cluster_num,
#       function(my_col) {

#         #Adding NA to the data when the a cluster is compared with itself
#         if (my_col==my_row) {
#           NA

#         #Checking if two clusters are not different enaugh
#         } else {
          
#           # Defining the mean density of each marker per cluster in named vectors
#           medianDens_cl1 <- sapply(markers,
#                                   function(marker) {
#                                     mean(DIN_matrix[clusters_object==my_row,markers])
#                                   }
#                                 )

#           medianDens_cl2 <- sapply(markers,
#                                   function(marker) {
#                                     mean(DIN_matrix[clusters_object==my_col,markers])
#                                   }
#                                 )

#           # Comparing values. Adding TRUE to the cluster_comp data frame if the clusters 
#           # are meaninfully different, else, adding FALSE.

#           # Check if a single marker is different
#           if (any(abs(medianDens_cl1 - medianDens_cl2) >= 0.003)) {
#             TRUE
          
#           # Check if there is an overal difference
#           } else if (sum(abs(medianDens_cl1 - medianDens_cl2)) >= 0.7 * 0.03 * length(markers) ) {
#             TRUE

#           # Returning FALSE if the conditions are not met
#           } else {
#             FALSE
#           }

#         }

#         # Returning the cluster_comp data frame
#         return(cluster_comp)

#       }
#     )
#   }
# }

# Running the cluster evaluation in parallel

#reviewed_optimal_clusters_list <- foreach(sample=names(optimal_clusters_list)) %dopar% {
  
  # Compare cell densities of the identified clusters
#  cluster_eval_df <- cluster_evaluation(DIN_matrix=DIN_matrices[[sample]], 
#                     clusters_object=optimal_clusters_list[[sample]], 
#                     markers=markers)

  # Merging clusters if required
#  for (cluster in rownames(cluster_eval_df)) {



#  }

#  list(sample, clusters)

#}

#
# SAVING OUTPUT FILES
#

# Generating output files
saveRDS(output_clusters_list, file="optimal_clusters.rds")