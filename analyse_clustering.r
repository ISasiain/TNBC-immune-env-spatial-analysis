rm(list = ls())
setwd("/home/isc/Spatial_immune_env")


#
# LOADING PACKAGES
#

#Loading packages for spatial analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")

#
# DEFINING VARIABLES TO RUN THE ANALYSIS
#

#Generate a list of the paths to the image objects
path_to_files <- "./data_from_suze/data/spiat/"
spe_objects <- list.files(path_to_files)

# Create a vector variable to store cell type(s) of interest
cell_type <- c("CD20", "CD3", "CD8")

# Generate a data frame to store information about clustering

#Generate colnames and rownames
column_names <- c("uid", "Clustered", "p_value", "Nr_of_clusters", "%_of_cells_in_clusters", "cluster_sizes", "free_cells")
row_names <- unname(
              sapply(spe_objects, 
                     function(tme) paste("BLOCK_", strsplit(tme, "_")[[1]][3], sep = "")
                     )
              )

#Generating empty data frame
clustering_metrics_df <- as.data.frame(matrix(NA, length(row_names), length(column_names)))
colnames(clustering_metrics_df) <- column_names
rownames(clustering_metrics_df) <- row_names


#
# DEFINING VARIABLES TO RUN THE ANALYSIS
#

for (tme in spe_objects) {
  
  tryCatch({
  
  my_tme <- readRDS(paste(path_to_files,tme,sep=""))
  uid <- paste("BLOCK_", strsplit(tme, "_")[[1]][3], sep="")
  
  
  # Analysing if the cell type of interest generates clusters
  check_clustering <- average_nearest_neighbor_index(my_tme, 
                                                        reference_celltypes = cell_type, 
                                                        feature_colname = "Cell.Type",
                                                        p_val = 01e-10)
  
  print(check_clustering)
  
  if (check_clustering$pattern == "Clustered") {
    
    #Identifying cell clusters
      clusters <- identify_neighborhoods(my_tme, 
                                       cell_types_of_interest = cell_type,
                                       feature_colname = "Cell.Type",
                                       radius =50,
                                       min_neighborhood_size = 35)

    if (length(na.omit(unique(clusters$Neighborhood)))==1) { #This means that only free cells have been identified
      
      
      clustering_metrics_df[uid,] <- c(
        
        "uid"=uid, #The sample identifier
        "Clustered"=FALSE, # The cells analysed do not generate clusters
        "p_value"=check_clustering$`p-value`, #Add the p_value measured
        "Nr_of_clusters"=NA, # This value is set to NA as no clusters have been identified
        "%_of_cells_in_clusters"=NA, # This value is set to NA as no clusters have been identified
        "cluster_sizes"=NA # This value is set to NA as no clusters have been identified
      )
      
    } else {
      
      #Calculate the composition of each cluster
      cluster_composition <- composition_of_neighborhoods(clusters, "Cell.Type")

      
      clustering_metrics_df[uid,] <- c(
        
        "uid"=uid, #The sample identifier
        "Clustered"=TRUE, # The cells analysed do not generate clusters
        "p_value"=check_clustering$`p-value`, #Add the p_value measured
        "Nr_of_clusters"=length(cluster_composition[,"Neighborhood"][cluster_composition[,"Neighborhood"]!="Free_cell"]), # 
        "%_of_cells_in_clusters"=1 - (cluster_composition[nrow(cluster_composition),"Number_of_cells"] / sum(cluster_composition[,"Number_of_cells"])),
        "cluster_sizes"=paste(head(cluster_composition$Number_of_cells, -1),collapse = ","),
        "free_cells"=cluster_composition$Number_of_cells[length(cluster_composition$Number_of_cells)] 
      )
      
    } 
    
  } else {
    
    clustering_metrics_df[uid,] <- c(
      
      "uid"=uid, #The sample identifier
      "Clustered"=FALSE, # The cells analysed do not generate clusters
      "p_value"=check_clustering$`p-value`, #Add the p_value measured
      "Nr_of_clusters"=NA, # This value is set to NA as no clusters have been identified
      "%_of_cells_in_clusters"=NA, # This value is set to NA as no clusters have been identified
      "cluster_sizes"=NA, # This value is set to NA as no clusters have been identified
      "free_cells"=NA # This value is set to NA as no clusters have been identified
    )
    
  }
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

saveRDS(clustering_metrics_df, file=paste0("clustering_metrics_", paste(cell_type, collapse = "-"), ".RData"))

