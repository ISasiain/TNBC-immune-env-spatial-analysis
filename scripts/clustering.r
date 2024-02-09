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

    make_option(c("-m", "--markers_for_clustering"), 
                type="character", 
                help="A comma separated list of the markers to be used to determine clusters should be entered here.",
                metavar="[PATH]"),

    make_option(c("-a", "--additional_clusters"), 
                type="character", 
                help="A comma separated list of the additional markers to analyse (not to define clusters) should be entered here.",
                metavar="[PATH]"),

    make_option(c("-M", "--minimum_distance_between_centroids"), 
                type="numeric",
                default=1*10^-4,
                help="The minimum mean distance between centroids of each sample to take the clustering into account.",
                metavar="[NUMBER]"),

    make_option(c("-D", "--density_threshold"), 
                type="numeric",
                default=3*10^-4,
                help="The minimum cell density in a sample to consider a sample feasible for the analysis.",
                metavar="[NUMBER]"),
    
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
original_samples <- names(DIN_matrices)

# Generating vector of markers 
markers <- strsplit(arguments$markers_for_clustering, ",")[[1]]
additional_markers <- strsplit(arguments$additional_clusters, ",")[[1]]


#
# CONFIGURATING PARALLELIZATION
#

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores)  
registerDoParallel(cl)  


#
# QC PREVIOUS TO CLUSTERING
#

# Filtering samples with low cell density
remove_samples <- unlist(foreach(sample=names(DIN_matrices)) %dopar% {
    
    #Determining total cell density of the markers for clustering
    density <- nrow(DIN_matrices[[sample]])/(pi*1500^2)

    #Determining the samples to be kept
    if (density < arguments$density_threshold) {
        return(sample)
    } else {
        return(NULL)
    }
    
})

# Print warning message about the dismissed samples
if (length(remove_samples) != 0) {

  cat("\n\n**** WARNING *****\n")
  cat("\nThe number of cells positive in any of the immune markers analysed of the following cores was considered too low, and therefore dismissed.\n")
  print(remove_samples)

}

saveRDS(remove_samples, "filtered_samples.rds")


# Removing samples with low cell density
DIN_matrices <- DIN_matrices[!(names(DIN_matrices) %in% remove_samples)]


#
# CLUSTERING BASED ON DIN VALUES
#

#Determining optimal number of clusters for each sample
optimal_number_of_clusters <- foreach(sample=names(DIN_matrices), .packages="NbClust") %dopar% {
  
  #If the density of some of the markers to be used for clustering is not defined
  #the sample will be ignored.
  if (all(markers %in% colnames(DIN_matrices[[sample]]))) {
    
    # Generate optimal number of clusters based on the specified markers
    tryCatch({(list(sample,
          length(unique(NbClust(
              DIN_matrices[[sample]][,markers],
              method="kmeans" #Using kmeans clustering
            )$Best.partition))
          )
        )},
           error=function(e) {} #Skip sample if error
        )
    

  }

}

#Rearraging list
num_clusters_list <- lapply(optimal_number_of_clusters, function(item) item[[2]])
names(num_clusters_list) <- sapply(optimal_number_of_clusters, function(item) item[[1]])

# Generating output files
saveRDS(num_clusters_list, file="number_optimal_clusters.rds")

#Remove NULL values from list
num_clusters_list <- num_clusters_list[!sapply(num_clusters_list, is.null)]



# Print warning message about the samples not used for clustering
not_clustered <- names(DIN_matrices)[!names(DIN_matrices) %in% names(num_clusters_list)]

if (length(not_clustered) != 0) {

  cat("\n\n**** WARNING *****\n")
  cat("\nThe following samples were not used for clustering as data for some of the markers was unavailable\n")
  print(not_clustered)

}



#Performing clustering
 output_clusters_list <- list()
 
 for (sample in names(num_clusters_list)) {
  
  #Perform kmeans clustering based on the optimal number of clusters
  tryCatch({

   output_clusters_list[[sample]] <- kmeans(DIN_matrices[[sample]][,markers], centers=num_clusters_list[[sample]])


  }, error=function(e){}) #Skip sample if error

 }


#
# DETERMINE TOO SIMILAR CLUSTERS FROM HOMOGENEOUS CORES
#

# Calculating the distance between centroids
dbc <- lapply(clusters, 
       function(cluster) {

                mean(dist( # Calculate distance between centers
                    cluster$centers # Getting centers of centroids
            ))
       }
    )

#Adding names to dbc list
names(dbc) <- names(clusters)

#Determining the samples with too similar clusters
homogeneous_samples <. foreach(sample=names(dbc)) {

    if (dbc[[sample]] < arguments$minimum_distance_between_centroids) {

        return(sample)

    } else {NULL}

}


# Print warining message about the dismissed clusters
cat("\n\nThe clusters detected in",
    sum(mean_dbc <= arguments$minimum_distance_between_centroids), 
    "samples were considered too similar, and therefore will be considered as homogeneous samples\n\n")


#
# CLUSTER ANNOTATION
#

# Generating list to store the output
annotated_samples <- vector("list", length(original_samples))
names(annotated_samples) <- original_samples


# ANNOTATING SAMPLES WITH LOW CELL DENSITY
annotated_samples[remove_samples] <- "Low cell density"


# ANNOTATING HOMOGENEOUS SAMPLES
annotated_samples[homogeneous_samples] <- foreach(sample=homogeneous_samples) %dopar% {

        # Generate matrix with p values obtained by Mann-Withney U test  
        p_value_greater <- matrix(NA, 
                                ncol=length(markers_for_clustering),
                                nrow=length(markers_for_clustering),
                                dimnames=list(markers_for_clustering,markers_for_clustering))


        # Filling the list of matrices
        for (row in markers_for_clustering) {
            for (col in markers_for_clustering) {

                if(row != col) {

                    p_value_greater[row, col] <- wilcox.test(
                        
                        DINs[[sample]][,row],
                        DINs[[sample]][,col],
                        exact=FALSE,
                        alternative = "greater")$p.value

                    }
                }

            }

  #Generating annotation
  list(

    type="Homogeneous",
    class = if (any(p_value_greater[row,] < 0.05)) {
      "TUMOUR"
     #} else if () {
     #  "IMMUNE"
     } else {
       "IMMUNE"
     }
   )
    
}


# ANNOTATING NON-HOMOGENEOUS SAMPLES
annotated_samples(original_samples[! original_samples %in% c(remove_samples, homogeneous_samples)]) <- foreach(sample=original_samples[! original_samples %in% c(remove_samples, homogeneous_samples)]) %dopar% {

    #Generating annotation
    list(
        type="Heterogeneous",
        clusters=lapply(sort(unique(unname(clusters[[sample]]$cluster))), function(clus) {

            #Checking if it is a tumour cluster
                        #Checking if it is a tumour cluster
            if (all(list_of_p_values[[sample]][[clus]]["p53", -which(colnames(list_of_p_values[[sample]][[clus]])=="p53")] < 1e-100)) {
                list(type = "TUMOUR"
                )

            #Checking if it is an immune cluster
           } else if (any(list_of_p_values[[sample]][[clus]][-which(rownames(list_of_p_values[[sample]][[clus]])=="p53"), "p53"] < 1e-100)) {

             # Determining abundance order of the markers in the immune and mixed clusters
                 order_vec <- rowSums(
                     list_of_p_values[[sample]][[clus]][
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53"), 
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53")] < 0.05,
                     na.rm = TRUE)

                 names(order_vec) <- rownames(list_of_p_values[[sample]][[clus]][
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53"), 
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53")])

                 groups <- sort(unique(order_vec), decreasing=TRUE)

                 for (element in 1:length(order_vec)) {

                   order_vec[element] <- which(groups == order_vec[element])

                 }



                list(
                    type = "IMMUNE",
                    order = order_vec
                )

            #If the previous conditions are not met, it will be considered as a mixed cluster 
            } else {
                # Determining abundance order of the markers in the immune and mixed clusters
                 order_vec <- rowSums(
                     list_of_p_values[[sample]][[clus]][
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53"), 
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53")] < 0.05,
                     na.rm = TRUE)

                 names(order_vec) <- rownames(list_of_p_values[[sample]][[clus]][
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53"), 
                         -which(colnames(list_of_p_values[[sample]][[clus]])=="p53")])

                 groups <- sort(unique(order_vec), decreasing=TRUE)

                 for (element in 1:length(order_vec)) {

                   order_vec[element] <- which(groups == order_vec[element])

                 }


                list(type = "MIXED"
                ,    order = order_vec
                )

            }

        })
    )

}


# Saving output list
saveRDS(annotated_samples, "annotated_samples.rds")