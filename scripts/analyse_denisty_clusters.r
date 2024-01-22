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

if (!require("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library(dplyr)


#
# PARSING COMMAND LINE ARGUMENTS
#


argument_list <- list(

    make_option(c("-d", "--DIN_matrix_list"), 
                type="character", 
                help="The path to the .rds file containing the list of DIN matrices of the analysed 
                        samples must be entered here.",
                metavar="[PATH]"),

    make_option(c("-l", "--clusters_list"), 
                type="character", 
                help="The path to the .rds file containing the list of the clustering output of all the analysed 
                        samples must be entered here.",
                metavar="[PATH]"),

    make_option(c("-m", "--minimum_distance_between_centroids"), 
                type="numeric",
                default=4.5*10^-4,
                help="The minimum mean distance between centroids of each sample to take the clustering into account.",
                metavar="[NUMBER]"),
    
    make_option(c("-a", "--markers_to_analyse"), 
              type="character",  
              help="Comma separated list of markers to analyse",
              metavar = "[COMMA_SEPARATAED_MARKERS]"))

arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program analysed the detected clusters based on DIN values."))


#
# READING INPUT FILES
#

#Annotation file
annotation <- read.csv(arguments$annotation, sep="\t")

#DIN list
DINs <- readRDS(arguments$DIN_matrix_list)

#Cluster list
clusters <- readRDS(arguments$clusters_list)

# Generating vector of markers to analyse
markers <- strsplit(arguments$markers_to_analyse, ",")[[1]]


#
# QUALITY CONTOL: DISMISSING TOO SIMILAR CLUSTERS
#

# Calculating the mean distance between centroids
mean_dbc <- sapply(clusters, 
       function(cluster) {

            mean( # Determining the mean
                dist( # Calculate distance between centers
                    cluster$centers # Getting centers of centroids
                )
            )
       }
    )

# Print warining message about the dismissed clusters
cat("\n\nThe clusters detected in",
    sum(mean_dbc <= arguments$minimum_distance_between_centroids), 
    "samples were considered too similar, and therefore will be dismissed\n")


# Eliminating dismissed clusters
clusters_list <- clusters_list[mean_dbc > arguments$minimum_distance_between_centroids]


#
# ESTABLISHING A SYSTEMATIC CLUSTER EVALUATION SYSTEM
#

# DATA FRAME TO STORE MEDIAN DENSITY DATA

# Generating data frame
cluster_dens_df <- data.frame()
colnames(cluster_dens_df) <- c("Cluster", "Sample", markers)


# Generating one row per each cluster

#Iterating through the samples
for(sample in names(clusters_list)) {

  #Iterating through the clusters
  for(cluster in 1:length(unique(clusters_list[[sample]]$cluster)))

  new_row <- c(
    cluster,
    sample,
    sapply(markers, function(marker) {
      median(
        DIN_matrices[[sample]][clusters_list[[sample]]==cluster,marker]
      )

  rbind(cluster_dens_df, new_row)

    })
  )

}

# DATA FRAME TO STORE CLUSTER CLASSIFICATION

# EXPLANATION OF CLUSTER CLASSIFICATION
# Very enriched (VE) <- 5th quintil of clusters (20% of the samples with highest value)
# Enriched (E) <- 4th quintil of clusters
# Neutral (N) <- 3rd quintil of clusters
# Diminished (D) <- 2nd quintil of clusters
# Very diminished (VD) <- 1st quintil of clusters

# Generate quintilles of the densities observed on the identified total clusters

#Generating data frame
classification_df <- data.frame()
colnames(classification_df) <- c("Cluster", "Sample", markers)

# Defining quintiles of each marker's density
quintiles_matrix <- matrix(NA, nrow=5, ncol=length(markers))
colnames(quintiles_matrix) <- markers
rownames(quintiles_matrix) <- c("0%", "20%", "40%", "60%", "80%", "100%")

# Defining quintiles
for (marker in markers) {

  quintiles_matrix[,marker] <- quantile(cluster_dens_df[,marker], probs = seq(0, 1, by = 0.2))

}

#Iterating through the rows of cluster_dens_df[
for(row in 1:nrow(cluster_dens_df)) {

  new_row <- c(
    cluster_dens_df[row, "Cluster"],
    scluster_dens_df[row, "Sample"],    
    
    sapply(markers, function(marker) {

      findInterval(cluster_dens_df[row, marker], quintiles_matrix[,marker])

    })
  )
}

# Mutate classification dataframes to replace interval numbers by factors.
value_mapping <- c("VD", "D", "N", "E", "VE")

classification_df <- classification %>% mutate_all(~factor(., levels = 1:5, labels = value_mapping))



#
# SAVING FILES
#

#Saving generated data
saveRDS(classification_df, file="cluster_classification.rds")
saveRDS(cluster_dens_df, file="cluster_median_density.rds")