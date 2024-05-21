#!/usr/bin/Rscript

# Description. This script generates cell clusters based on DIN values.
# Usage. clustering.r -d [path] -m [comma separated list] -t [name] -l [number] -M [number] -D [number] -c [number] -n [name] -p [path]
# Arguments: -d [path] The path to the .rds file containing the list of DIN matrices of the analysed samples must be entered here.
#            -m [comma separated list] A comma separated list of the markers to be used to determine clusters should be entered here.
#            -a [comma separated list] A comma separated list of the additional markers to analyse (not to define clusters) should be entered here.
#            -l [number] The maximum number of expected clusters must be entered here.
#            -t [name] The name of the marker to use as a tumour markers must be entered here.
#            -M [number] The minimum mean distance between centroids of each sample to take the clustering into account.
#            -D [number] The minimum cell density in a sample to consider a sample feasible for the analysis.
#            -c [number] Number of cores to be used to run the program [default 1]
#            -n [name] Name of the output files [default output_clustering]
#            -p [path] Path of the output files [default ./]
# Output: Annotated samples and clustered cells in the output directory.


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

    make_option(c("-a", "--additional_markers"), 
                type="character", 
                default=NA,
                help="A comma separated list of the additional markers to analyse (not to define clusters) should be entered here.",
                metavar="[PATH]"),

    make_option(c("-l", "--max_number_of_clusters"), 
                type="numeric", 
                help="The maximum number of expected clusters must be entered here.",
                default=5,
                metavar="[PATH]"),

    make_option(c("-t", "--tumour_marker"), 
                type="character", 
                help="The name of the marker to use as a tumour markers must be entered here.",
                metavar="[PATH]"),

    make_option(c("-M", "--minimum_distance_between_centroids"), 
                type="numeric",
                default=1.20,
                help="The minimum mean distance between centroids of each sample to take the clustering into account.",
                metavar="[NUMBER]"),

    make_option(c("-D", "--density_threshold"), 
                type="numeric",
                default=2.8*10^-4,
                help="The minimum cell density in a sample to consider a sample feasible for the analysis.",
                metavar="[NUMBER]"),
    
    make_option(c("-c", "--cores"), 
              type="integer", 
              default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[NUMBER]"),

    make_option(c("-n", "--output_name"), 
              type="character", 
              default="output_clustering",  
              help="Name of the output files [default %default]",
              metavar = "[NAME]"),

    make_option(c("-p", "--path_output"), 
              type="character", 
              default="./",  
              help="Path of the output files [default %default]",
              metavar = "[PATH]")
              
              )

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
additional_markers <- strsplit(arguments$additional_markers, ",")[[1]]


#
# CONFIGURATING PARALLELIZATION
#

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores, outfile="clustering.log")  
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

# Checking and removing samples that do not contain any tumour marker positive cell
without_tumour_samples <- c()
for (sample in names(DIN_matrices)) {
    
    if (!arguments$tumour_marker %in% DIN_matrices[[sample]]$Phenotype) {
        without_tumour_samples <- c(without_tumour_samples, sample)
    }
    
}

# Print warning message about the dismissed samples
if (length(without_tumour_samples) != 0) {

  cat("\n\n**** WARNING *****\n")
  cat("\nThe following cores did not contain any tumour cell; therefore will be dismissed.\n")
  print(without_tumour_samples)

}


# Removing samples with low cell density or without tumour cells
DIN_matrices <- DIN_matrices[!(names(DIN_matrices) %in% remove_samples)]
DIN_matrices <- DIN_matrices[!(names(DIN_matrices) %in% without_tumour_samples)]

#
# CLUSTERING
#

cat("\nFinding optimal number of clusters...\n")


#Determining optimal number of clusters for each sample
optimal_number_of_clusters <- foreach(sample=names(DIN_matrices), .packages="NbClust") %dopar% {

    if (arguments$max_number_of_clusters > 2) {

        # Generate optimal number of clusters based on the specified markers that are included in the sample
        tryCatch({(list(sample,
            length(unique(NbClust(
                DIN_matrices[[sample]][,markers[markers %in% colnames(DIN_matrices[[sample]])]],
                max.nc=arguments$max_number_of_clusters, #Maximum number of clusters
                method="kmeans" #Using kmeans clustering
            )$Best.partition))
            )
        )},
            error=function(e) {} #Skip sample if error
        )
    } else {list(sample, 2)} # If the maximum number of clusters is two there is no need to estimate the optimal number of clusters
}

#Rearraging list
num_clusters_list <- lapply(optimal_number_of_clusters, function(item) item[[2]])
names(num_clusters_list) <- sapply(optimal_number_of_clusters, function(item) item[[1]])

#Remove NULL values from list
num_clusters_list <- num_clusters_list[!sapply(num_clusters_list, is.null)]



#Performing clustering
 output_clusters_list <- list()
 
 for (sample in names(num_clusters_list)) {
  
  #Perform kmeans clustering based on the optimal number of clusters
  tryCatch({

   output_clusters_list[[sample]] <- kmeans(DIN_matrices[[sample]][,markers[markers %in% colnames(DIN_matrices[[sample]])]], centers=num_clusters_list[[sample]])


  }, error=function(e){}) #Skip sample if error

 }


#
# DETERMINE TOO SIMILAR CLUSTERS FROM HOMOGENEOUS CORES
#

# Calculating the distance between centroids
dbc <- lapply(output_clusters_list, 
       function(cluster) {

                mean(dist( # Calculate distance between centers
                    cluster$centers, # Getting centers of centroids
                    method="canberra" # Using canberra distance
            ))
       }
    )

#Adding names to dbc list
names(dbc) <- names(output_clusters_list)


#Determining the samples with too similar clusters
homogeneous_samples <- names(dbc)[dbc <= arguments$minimum_distance_between_centroids]

# Print warining message about the dismissed clusters
cat("\n\nThe clusters detected in",
    sum(dbc <= arguments$minimum_distance_between_centroids), 
    "samples were considered too similar, and therefore will be considered as homogeneous samples\n\n")


#
# CLUSTER ANNOTATION
#

# Adding additional markers to the list of markers
if (length(na.omit(additional_markers)) != 0) {
    markers <- c(markers, additional_markers)
}


# Generating list to store the output
annotated_samples <- vector("list", length(original_samples))
names(annotated_samples) <- original_samples


# ANNOTATING SAMPLES WITH LOW CELL DENSITY
annotated_samples[remove_samples] <- "Low cell density"

# ANNOTATING SAMPLES WITHOUT TUMOUR CELLS
annotated_samples[without_tumour_samples] <- "Without tumour cells"

# ANNOTATING HOMOGENEOUS SAMPLES
for (sample in homogeneous_samples) {
        
        # Generate matrix with p values obtained by Mann-Withney U test  
        p_value_greater <- matrix(NA, 
                                ncol=length(markers[markers %in% colnames(DIN_matrices[[sample]])]),
                                nrow=length(markers[markers %in% colnames(DIN_matrices[[sample]])]),
                                dimnames=list(markers[markers %in% colnames(DIN_matrices[[sample]])],markers[markers %in% colnames(DIN_matrices[[sample]])]))


        # Filling the list of matrices
        for (row in markers[markers %in% colnames(DIN_matrices[[sample]])]) {
            for (col in markers[markers %in% colnames(DIN_matrices[[sample]])]) {

                if(row != col) {

                    p_value_greater[row, col] <- wilcox.test(
                        
                        DIN_matrices[[sample]][,row],
                        DIN_matrices[[sample]][,col],
                        exact=FALSE,
                        alternative = "greater")$p.value

                                        }
                                    }

                                }


                      #Generating annotation
                      annotated_samples[[sample]] <- list(

                        type="Homogeneous",
                        class = if (all(p_value_greater[arguments$tumour_marker, -which(colnames(p_value_greater)==arguments$tumour_marker)] < 1e-100)) {
                                    list(type = "TUMOUR"
                                    )

                                #Checking if it is an immune cluster
                               } else if (any(p_value_greater[-which(rownames(p_value_greater)==arguments$tumour_marker), arguments$tumour_marker] < 1e-100)) {

                                  #Determining abundance order of the markers[markers %in% colnames(DIN_matrices[[sample]])] in the immune and mixed clusters
                                     order_vec <- colSums(
                                         p_value_greater[
                                             -which(colnames(p_value_greater)==arguments$tumour_marker), 
                                             -which(colnames(p_value_greater)==arguments$tumour_marker)] < 0.05,
                                         na.rm = TRUE)

                                     names(order_vec) <- rownames(p_value_greater[
                                             -which(colnames(p_value_greater)==arguments$tumour_marker), 
                                             -which(colnames(p_value_greater)==arguments$tumour_marker)])

                                     groups <- sort(unique(order_vec), decreasing=TRUE)

                                     for (element in 1:length(order_vec)) {

                                       order_vec[element] <- which(groups == order_vec[element])

                                     }



                                    list(
                                        type = "IMMUNE"
                                        ,order = order_vec
                                    )

                                #If the previous conditions are not met, it will be considered as a mixed cluster 
                                } else {
                                    # Determining abundance order of the markers in the immune and mixed clusters
                                     order_vec <- rowSums(
                                         p_value_greater[
                                             -which(colnames(p_value_greater)==arguments$tumour_marker), 
                                             -which(colnames(p_value_greater)==arguments$tumour_marker)] < 0.05,
                                         na.rm = TRUE)

                                     names(order_vec) <- rownames(p_value_greater[
                                             -which(colnames(p_value_greater)==arguments$tumour_marker), 
                                             -which(colnames(p_value_greater)==arguments$tumour_marker)])

                                     groups <- sort(unique(order_vec), decreasing=TRUE)

                                     for (element in 1:length(order_vec)) {

                                       order_vec[element] <- which(groups == order_vec[element])

                                     }


                                    list(type = "MIXED"
                                        ,order = order_vec
                                    )
                          }
                        )
                        
                    }


                    # ANNOTATING NON-HOMOGENEOUS SAMPLES

list_of_p_values <- list()

for (sample in original_samples[! original_samples %in% c(remove_samples, homogeneous_samples)]) {

  # Generating list of matrices with p values obtained by Mann-Withney U test

    list_of_p_values <- 

          lapply(sort(unique(unname(output_clusters_list[[sample]]$cluster))), function(clus) {
          
          # Generate list of matrices with p values obtained by Mann-Withney U test and effect size (difference of median)
          p_value_greater <- matrix(NA, 
                                  ncol=length(markers[markers %in% colnames(DIN_matrices[[sample]])]),
                                  nrow=length(markers[markers %in% colnames(DIN_matrices[[sample]])]),
                                  dimnames=list(markers[markers %in% colnames(DIN_matrices[[sample]])],markers[markers %in% colnames(DIN_matrices[[sample]])]))
          

          # Filling the list of matrices

          for (row in markers[markers %in% colnames(DIN_matrices[[sample]])]) {
              for (col in markers[markers %in% colnames(DIN_matrices[[sample]])]) {
                  if(row != col) {
                      p_value_greater[row, col] <- wilcox.test(
                          
                          DIN_matrices[[sample]][output_clusters_list[[sample]]$cluster==clus,row],
                          DIN_matrices[[sample]][output_clusters_list[[sample]]$cluster==clus,col],
                          exact=FALSE,
                          alternative = "greater")$p.value

                      }
                  }
              }
              p_value_greater
          })


    #Generating annotation
    annotated_samples[[sample]] <- list(
        type="Heterogeneous",
        clusters=lapply(sort(unique(unname(output_clusters_list[[sample]]$cluster))), function(clus) {

            #Checking if it is a tumour cluster
            if (all(list_of_p_values[[clus]][arguments$tumour_marker, -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker)] < 1e-100)) {
                list(type = "TUMOUR")

            #Checking if it is an immune cluster
            } else if (any(list_of_p_values[[clus]][-which(rownames(list_of_p_values[[clus]])==arguments$tumour_marker), arguments$tumour_marker] < 1e-100)) {

             # Determining abundance order of the markers in the immune and mixed clusters
                 order_vec <- rowSums(
                     list_of_p_values[[clus]][
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker), 
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker)] < 0.05,
                     na.rm = TRUE)

                 names(order_vec) <- rownames(list_of_p_values[[clus]][
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker), 
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker)])

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
                     list_of_p_values[[clus]][
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker), 
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker)] < 0.05,
                     na.rm = TRUE)

                 names(order_vec) <- rownames(list_of_p_values[[clus]][
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker), 
                         -which(colnames(list_of_p_values[[clus]])==arguments$tumour_marker)])

                 groups <- sort(unique(order_vec), decreasing=TRUE)

                 for (element in 1:length(order_vec)) {

                   order_vec[element] <- which(groups == order_vec[element])

                 }


                list(type = "MIXED",
                    order = order_vec
                )

            }

        })
    )

}


# Saving output list
saveRDS(annotated_samples, paste0(arguments$path_output, arguments$output_name, ".rds"))
saveRDS(output_clusters_list,  paste0(arguments$path_output, arguments$output_name, "_clustered_cells.rds"))


cat("\n\n****************\n")
cat("Process finished!\n")
cat("****************\n\n")    

# Stop parallelization
stopCluster(cl)