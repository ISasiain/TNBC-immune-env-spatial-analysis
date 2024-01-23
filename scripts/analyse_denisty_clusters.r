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

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)


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
                default=4*10^-4,
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
#annotation <- read.csv(arguments$annotation, sep="\t")

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

# Generate histogram with the dbc of all the clusters detected

ggplot(data=data.frame(x=mean_dbc), aes(x=x)) +
  geom_histogram() +
  geom_vline(xintercept=arguments$minimum_distance_between_centroids, linetype="dashed", color="darkolivegreen3", size=2) +
  xlab("Mean distance between clusters") +
  ylab("Number of samples") +
  theme_classic()

ggsave("mean_dbc_of_samples.pdf", device="pdf")


# Print warining message about the dismissed clusters
cat("\n\nThe clusters detected in",
    sum(mean_dbc <= arguments$minimum_distance_between_centroids), 
    "samples were considered too similar, and therefore will be dismissed\n\n")


# Eliminating dismissed clusters
clusters <- clusters[mean_dbc > arguments$minimum_distance_between_centroids]

#
# ESTABLISHING A SYSTEMATIC CLUSTER EVALUATION SYSTEM
#

# DATA FRAME TO STORE MEDIAN DENSITY DATA

# Generating data frame
# Create an empty data frame with specified column names
cluster_dens_df <- data.frame()

# Generating one row per each cluster

# Iterating through the samples
for (sample in names(clusters)) {

  # Iterating through the clusters
  for (cluster in 1:length(unique(clusters[[sample]]$"cluster"))) {

    # Creating a new row
    new_row <- c(
      cluster,
      sample,
      sapply(markers, function(marker) {
        
        #Using tryCatch to prevent the script from failing if one of the markers is not included 
        #in one of the DIN matrices
        tryCatch({
          median(
            DINs[[sample]][clusters[[sample]]$"cluster" == cluster, marker]
          )
        },error=function(e) {NaN})

      })
    )

    # Appending the new row to the data frame
    cluster_dens_df <- rbind(cluster_dens_df, new_row)
  }
}

#Setting colnames of the cluster_dens_df data frame
colnames(cluster_dens_df) <- c("Cluster", "Sample", markers)


# Generate histogram with the cel densities of the clusters detected
for (marker in markers) {

#Change data type of densities
to_plot <- na.omit(as.numeric(cluster_dens_df[,marker]))

#pdf(paste0(marker,"_density_of_clusters.pdf"))

ggplot(data=data.frame(x=to_plot), aes(x=x)) +
  geom_histogram() +
  xlab("Median density in the cluster") +
  ylab("Number of samples") +
  theme_classic()

#hist(cluster_dens_df[,marker])
ggsave(paste0(marker,"_density_of_clusters.pdf"), device="pdf")

}

# Plot scatterplots of densities of each marker against all the others
for (marker1 in markers) {

  for (marker2 in markers) {
    
    to_plot <- na.omit(data.frame(
    as.numeric(cluster_dens_df[,marker1]),
    as.numeric(cluster_dens_df[,marker2])
    ))

    colnames(to_plot) <- c("marker1_DIN", "marker2_DIN")



    ggplot(data=to_plot, aes(x=marker1_DIN, y=marker2_DIN)) +
      geom_point() +
      xlab(paste0("Median DIN of ", marker1)) +
      ylab(paste0("Median DIN of ", marker2)) +
      theme_classic()

    ggsave(paste0(marker1,"_vs_", marker2, "_clusters_scatterplot.pdf"), device="pdf")

  }
}


# DATA FRAME TO STORE CLUSTER CLASSIFICATION

# EXPLANATION OF CLUSTER CLASSIFICATION
# Very enriched (VE) <- 5th quintil of clusters (20% of the samples with highest value)
# Enriched (E) <- 4th quintil of clusters
# Neutral (N) <- 3rd quintil of clusters
# Diminished (D) <- 2nd quintil of clusters
# Very diminished (VD) <- 1st quintil of clusters

# Generate quintilles of the densities observed on the identified total clusters

# Generating an empty data frame with specified column names
classification_df <- data.frame()

# Defining quintiles of each marker's density
quintiles_matrix <- matrix(, nrow = 6, ncol = length(markers))
colnames(quintiles_matrix) <- markers
rownames(quintiles_matrix) <- c("0%", "20%", "40%", "60%", "80%", "100%")

# Defining quintiles
for (marker in markers) {
  quintiles_matrix[, marker] <- quantile(as.numeric(cluster_dens_df[, marker]), probs = seq(0, 1, by = 0.2), na.rm=TRUE)
}


# Iterating through the rows of cluster_dens_df
for (row in 1:nrow(cluster_dens_df)) {
  new_row <- c(
    cluster_dens_df[row, "Cluster"],
    cluster_dens_df[row, "Sample"],
    sapply(markers, function(marker) {
      findInterval(cluster_dens_df[row, marker], quintiles_matrix[, marker])
    })
  )

  # Adding the new row to the classification_df data frame
  classification_df <- rbind(classification_df, new_row)
}

#Setting colnames of the data frame
colnames(classification_df) <- c("Cluster", "Sample", markers)

# Mutate classification_df to replace interval numbers by factors.
value_mapping <- c("VD", "D", "N", "E", "VE")

classification_df <- classification_df %>%
  mutate_at(vars(-Cluster, -Sample), ~factor(., levels = 1:5, labels = value_mapping))



#
# SAVING FILES
#

#Saving generated data
saveRDS(classification_df, file="cluster_classification.rds")
saveRDS(cluster_dens_df, file="cluster_median_density.rds")