rm(list = ls())
setwd("/home/isc/Spatial_immune_env/Spatial_immune_env_analysis")

#Loading packages for spatial analysis
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

#Uploading some example images to analyse
J10_B1 <- readRDS("../data_from_suze/data/spiat/PD31171a_BLOCK_1|J10_spiat.rds")
E3_B2 <- readRDS("../data_from_suze/data/spiat/PD31042a_BLOCK_3|I5_spiat.rds")
B7_B3 <- readRDS("../data_from_suze/data/spiat/PD36039a_BLOCK_3|B7_spiat.rds")
L5_B4 <- readRDS("../data_from_suze/data/spiat/PD36041a_BLOCK_4|L5_spiat.rds")
E7_B5 <- readRDS("../data_from_suze/data/spiat/PD31179a_BLOCK_5|E7_spiat.rds")
B2_B5 <- readRDS("../data_from_suze/data/spiat/PD31111a_BLOCK_5|B2_spiat.rds")
F4_B2 <- readRDS("../data_from_suze/data/spiat/PD35967a_BLOCK_2|F4_spiat.rds")
J7_B3 <- readRDS("../data_from_suze/data/spiat/PD35982a_BLOCK_3|J7_spiat.rds")
A2_B4 <- readRDS("../data_from_suze/data/spiat/PD31029a_BLOCK_4|A2_spiat.rds")
K10_B3 <- readRDS("../data_from_suze/data/spiat/PD31030a_BLOCK_3|K10_spiat.rds")
F10_B1 <- readRDS("../data_from_suze/data/spiat/PD36046a_BLOCK_1|F10_spiat.rds")
#Plotting cell types
plot_cell_categories(J10_B1, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(E3_B2, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(B7_B3, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(L5_B4, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(E7_B5, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(B2_B5, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(F4_B2, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(J7_B3, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(F4_B2, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(A2_B4, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(K10_B3, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))
plot_cell_categories(F10_B1, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))


unlist(table(J10_B1$Cell.Type) / (2*pi*1500))["CD20"]
#
# DEFINING FUNCTIONS TO DETERMINE PARAMETERS
#

#Analyising cells in neghborhood
neighbour_cells <- function(spe.object, radius, feature.colname="Cell.Type") {
  
  #Define cell types in the spe object
  cell_types <- unique(spe.object$Cell.Type)
  
  #Create data frame to store the data
  neighbours_df <- matrix(NA, length(cell_types), length(cell_types), dimnames = list(cell_types, cell_types))
  
  #Fill the dataframe
  
  for(reference in cell_types) {
    
    neighbours_df[reference,] <- sapply(cell_types, function(my_cell_type) (average_percentage_of_cells_within_radius(spe_object = spe.object, feature_colname = feature.colname, reference_celltype = reference, target_celltype = my_cell_type, radius = radius)))
    
  }
  
  #The rows are reference cells, and columns target cells
  #The sum of the row must be 100
  return(neighbours_df)
  
}

#Creating a function to check the neighborhood proportions minus the expected proportion
neighbour_props_vs_expected <- function(spe.object, rad, feature.colname="Cell.Type") {
  
  #Define cell types in the spe object
  cell_types <- unique(spe.object$Cell.Type)
  
  #Generating dataframe with neighbour cells
  neighbour_props <- neighbour_cells(spe.object, radius = rad, feature.colname = feature.colname)

  
  #Substracting the expected proportion
  for (col in cell_types) {
    neighbour_props[,col] <- neighbour_props[,col] - 100*(table(spe.object$Phenotype)[col]/length(spe.object$Phenotype))
  }
  
  return(neighbour_props)
  
}


#Create function to plot evolution of neighborhood proportions minus the expected proportion per cell type
plot_evolution_props_vs_expected <- function(spe, min_radius, max_radius, step_radius) {
  
  # Define cell types in the spe object
  cell_types <- unique(spe$Cell.Type)
  
  #Generate objects to store values Using a 3D array
  first_dim_reference <- cell_types
  second_dim_target <- cell_types
  third_dim_radius <- seq(min_radius, max_radius, step_radius)
  
  metrics <- array(data = NA,
                   dim = c(length(first_dim_reference), length(second_dim_target), length(third_dim_radius)),
                   dimnames = list(first_dim_reference,
                                   second_dim_target,
                                   as.character(third_dim_radius))
  )

  # Iterate through radii values
  for (r in seq(min_radius, max_radius, step_radius)) {
  
    metrics_for_radius <- neighbour_props_vs_expected(spe.object = spe, rad = r)
    
    for (cell_ref in cell_types) {
      for (cell_tar in cell_types)  {
        metrics[cell_ref, cell_tar, as.character(r)] <- metrics_for_radius[cell_ref, cell_tar]
      }
    }
  }

  #Generate list to store plots
  plot_list <- list()
  
  #Generate a plot per each reference cell type
  for (ref_cell in cell_types) {
    #Obtaining data from 3D array
    to_plot <- as.data.frame(metrics[ref_cell, , ])
    to_plot$"target_cell" <- rownames(to_plot)

    #Reshaping dataframe to be plotted using ggplot
    to_plot <- to_plot %>% gather(key="radius", value = "difference", -"target_cell")
    
    #Plotting a graph per ref cell
    plot_list [[ref_cell]] <- ggplot(to_plot, aes(x = as.numeric(radius), y = difference, group = as.factor(target_cell), color = as.factor(target_cell))) +
      geom_line(linewidth=3) +
      geom_hline(yintercept=0) +
      theme_minimal() + 
      labs(x = "Radius", y = "Δ CIN", title = paste("Reference Cell Phenotype:", ref_cell, sep=" ")) +
      guides(color = guide_legend(title = "Target Cell Phenotype"))
  }
  
  return(plot_list)
  
}


#
# DEFINE FUNCTION TO PLOT COLOCALIZATION OF CELLS
#


plot_colocalization <- function(spe_object, radius, feature_colname) {
  
  #Defining cell phenotypes
  cell_types <- unique(spe_object$Cell.Type)
  
  #Get cell proportions in the whole sample
  cell_props <- 100*(table(spe_object$Phenotype)/length(spe_object$Phenotype))
  
  #Determine observed - expected CIN
  deltaCIN <- neighbour_props_vs_expected(spe.object = spe_object, rad = radius, feature.colname = feature_colname)
  
  #Generating plot. IMPROVE APPEARANCE!!
  plotmat(deltaCIN,
          lwd = 3,
          box.lwd = 2,
          box.size = 0.1,
          box.type = "circle",
          box.col = c("skyblue", "lightgreen", "salmon", "orange"),
          arr.length = 0.5,
          arr.width = 0.2,
          cex.txt = 1,
          box.prop = 0.5,
          main = "Colocalization plot")
}

#
# DEFINE FUNCTION TO CALCULATE CELL DENSITY
#

cell_density <- function(spe_object, radius_length) {
  
  # Counting the number of cells in the sample
  nCells <- unlist(table(spe_object$Phenotype))

  # Divide the number of cells by the area of the sample
  return(nCells/(pi*radius_length^2))
  
}

cell_density(J10_B1, 1500)


cell_density_in_neighborhood <- function(spe_object, reference_cell, target_cell, radius) {
  
  # Combine cell ids, phenotypes, and coordinates of cells in a single data frame
  spatial_data <- cbind(spe_object@colData[, c("Cell.ID", "Phenotype")], spatialCoords(spe_object))
  
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

  

spatialCoords(J10_B1)
J10_B1@int_colData
cbind(J10_B1@colData[,c("Cell.ID", "Phenotype")], spatialCoords(J10_B1))

#  
# CALCULATING PARAMETERS
#

# Average pairwise distance (APD)

#Calculating distances between cells and plotting them
APD <- calculate_pairwise_distances_between_celltypes(J10_B1, cell_types_of_interest = c("p53"), feature_colname = "Cell.Type")
plot_cell_distances_violin(APD)

#Determine and plot summary statistics
APD_sum <- calculate_summary_distances_between_celltypes(APD)
plot_distance_heatmap(APD_sum)

# Average minimum distance (AMD)

#Calculating distances between cells and plotting them
AMD <-calculate_minimum_distances_between_celltypes(J10_B1, cell_types_of_interest = c("CD3", "p53"), feature_colname = "Cell.Type")
plot_cell_distances_violin(AMD)

#Determine and plot summary statistics
AMD_sum <- calculate_summary_distances_between_celltypes(AMD)
plot_distance_heatmap(AMD_sum)

# Cells in the neighbourhood
neighbour_cells(L5_B4, radius = 100)

# Percentage of cells in a radius vs radius (ADAPT YLIM DEPENDING ON THE CELLS COMPARED)
for (r in c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400)) {
  if (r==10) {
    plot(r, average_percentage_of_cells_within_radius(L5_B4, reference_celltype = "CD20", target_celltype = "CD3", feature_colname = "Cell.Type", radius = r), xlim=c(0,400), ylim=c(0,50))
  } else {
    points(r, average_percentage_of_cells_within_radius(L5_B4, reference_celltype = "CD20", target_celltype = "CD3", feature_colname = "Cell.Type", radius = r))
    }
  }

a <- neighbour_props_vs_expected(B2_B5, rad = 100)


#Location of cells by radii resolved CIN - Expected prop
a <- plot_evolution_props_vs_expected(B2_B5, 100, 2000, 50)
par(mfrow=c(1,4))
a$CD3
a$CD20
a$CD8
a$p53


# Mixing score and normalized mixing score  (MS and NMS)
for (r in c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1500, 2000)) {
 
  print(mixing_score_summary(F4_B2, 
                       reference_celltype = "CD3", 
                       target_celltype = "CD20", 
                       feature_colname = "Cell.Type",
                       radius = r)$Normalised_mixing_score)
}


#K cross function analysis
df_cross <- calculate_cross_functions(J10_B1, 
                                      method = "Kcross", 
                                      cell_types_of_interest = c("p53", "CD8"), 
                                      feature_colname = "Cell.Type", 
                                      dist = 200)
AUC_of_cross_function(df_cross) #Only works if the functions do not cross

crossing_of_crossK(df.cross = df_cross) #Only works if the functions cross

# Entropy

general_entropy <- calculate_entropy(F4_B2, 
                                     cell_types_of_interest = immune_cells, 
                                     feature_colname = "Cell.Type",
                                     radius = 1000)
print(general_entropy)

gradient_positions <- seq(100, 5000, 50)
gradient_entropy <- 
  compute_gradient(F4_B2, radii = gradient_positions, 
                   FUN = calculate_entropy,  cell_types_of_interest = c("CD8","p53"),
                   feature_colname = "Cell.Type")
length(gradient_entropy)
head(gradient_entropy)


gradient_results <- entropy_gradient_aggregated(F4_B2, radii = gradient_positions, 
                                                cell_types_of_interest = c("CD3","p53"),
                                                feature_colname = "Cell.Type")
plot(gradient_positions,gradient_results$gradient_df[1, 3:ncol(gradient_results$gradient_df)])


#Entropy per zone
entropy_grid <- grid_metrics(B2_B5,
                         FUN = calculate_entropy,
                         n_split = 40,
                         feature_colname="Cell.Type",
                         cell_types_of_interest = c("CD3", "CD20"))

#Identifying tumor cell sections
borders <- identify_bordering_cells(F4_B2, reference_cell = "p53", feature_colname = "Cell.Type")
distances <- calculate_distance_to_margin(borders)
immune_cells <- c("CD3", "CD8", "CD20")

struc <- define_structure(
  distances, cell_types_of_interest = immune_cells, 
  feature_colname = "Cell.Type", n_margin_layers = 4)

plot_cell_categories(struc, feature_colname = "Structure")

immune_props <- calculate_proportions_of_cells_in_structure(structure, immune_cells, feature_colname = "Cell.Type")



calculate_minimum_distances_between_celltypes(B2_B5, feature_colname = "Cell.Type", cell_types_of_interest = c("p53"))


#Identifying cell neighborhoods
clusters <- identify_neighborhoods(A2_B4,
                                   method = "hierarchical",
                                   cell_types_of_interest = c("p53"),
                                   feature_colname = "Cell.Type",
                                   radius = 50,
                                   min_neighborhood_size = 65)

cluster_composition <- composition_of_neighborhoods(clusters, "Cell.Type")

number_of_clusters <- length(cluster_composition[,"Neighborhood"][cluster_composition[,"Neighborhood"]!="Free_cell"])
cells_in_clusters <- 1 - (cluster_composition[nrow(cluster_composition),"Number_of_cells"] / sum(cluster_composition[,"Number_of_cells"]))


average_nearest_neighbor_index(A2_B4, 
                               reference_celltypes = "p53", 
                               feature_colname = "Cell.Type",
                               p_val = 01e-10)



function (spe_object, method = "hierarchical", cell_types_of_interest, 
          radius, min_neighborhood_size = 10, k = 100, feature_colname, 
          no_pheno = NULL) 
{
  Cell.X.Position <- Cell.Y.Position <- Cluster <- Xpos <- Ypos <- NULL
  plot_clusters <- TRUE
  formatted_data <- get_colData(spe_object)
  if (!is.null(no_pheno)) {
    formatted_data <- formatted_data[formatted_data[, feature_colname] != 
                                       no_pheno, ]
  }
  if (!is.null(cell_types_of_interest)) {
    formatted_data <- formatted_data[formatted_data[, feature_colname] %in% 
                                       cell_types_of_interest, ]
  }
  if (nrow(formatted_data) == 0) {
    message("There are no cells in data/no cells for the phenotypes of interest")
    formatted_data <- get_colData(spe_object)
    formatted_data$Cluster <- NA
    plot_clusters <- FALSE
  }
  else {
    if (method == "hierarchical") {
      rownames(formatted_data) <- formatted_data$Cell.ID
      sim_close <- -apcluster::negDistMat(formatted_data[, 
                                                         c("Cell.X.Position", "Cell.Y.Position")])
      sim_close[sim_close > radius] <- NA
      sim_close[sim_close == 0] <- NA
      sim_close <- sim_close[apply(sim_close, 1, function(x) {
        if (sum(is.na(x)) == length(x)) {
          return(FALSE)
        }
        else {
          return(TRUE)
        }
      }), ]
      if (is.null(dim(sim_close))) {
        formatted_data$Cluster <- NA
        plot_clusters <- FALSE
      }
      else {
        if (nrow(sim_close) == 0) {
          formatted_data[formatted_data[[feature_colname]] %in% 
                           cell_types_of_interest, "Cluster"] <- "Free_cell"
          plot_clusters <- FALSE
          message("There are no clusters detected in this image. All cells of interest are free cells.")
        }
        else {
          sim_close <- sim_close[, apply(sim_close, 
                                         2, function(x) {
                                           if (sum(is.na(x)) == length(x)) {
                                             return(FALSE)
                                           }
                                           else {
                                             return(TRUE)
                                           }
                                         })]
          cells_in_cluster <- rownames(sim_close)
          sim_close <- ifelse(is.na(sim_close), 1, 0)
          print(sim_close)
          if (nrow(sim_close) != 0 & ncol(sim_close) != 
              0) {
            h <- stats::hclust(stats::as.dist(sim_close), 
                               method = "single")
            local_clusters <- stats::cutree(h, h = 0.5)
            formatted_data$Cluster <- as.character(local_clusters[match(formatted_data$Cell.ID, 
                                                                        names(local_clusters))])
          }
          else {
            formatted_data$Cluster <- NA
          }
          summarised_data <- formatted_data %>% group_by(Cluster) %>% 
            summarise(n = n())
          big_clusters <- summarised_data[summarised_data$n > 
                                            min_neighborhood_size, "Cluster"]$Cluster
          print(big_clusters)
          if (length(big_clusters) == 0) {
            formatted_data[formatted_data[[feature_colname]] %in% 
                             cell_types_of_interest, "Cluster"] <- "Free_cell"
            plot_clusters <- FALSE
            message("There are no clusters detected in this image. All cells of interest are free cells.")
          }
          else {
            formatted_data[formatted_data$Cluster %in% 
                             big_clusters, "size"] <- "larger"
            cluster_ids <- unique(formatted_data[formatted_data$size == 
                                                   "larger", "Cluster"])
            n_cluster <- length(cluster_ids)
            n <- 1
            for (cluster_id in cluster_ids) {
              if (!is.na(cluster_id)) {
                formatted_data[which(formatted_data$Cluster == 
                                       cluster_id), "new_cluster"] <- n
                n <- n + 1
              }
            }
            formatted_data$Cluster <- formatted_data$new_cluster
            formatted_data$Cluster <- as.character(formatted_data$Cluster)
            formatted_data$new_cluster <- NULL
            formatted_data$size <- NULL
          }
        }
      }
    }
    else if (method == "dbscan") {
      cell_cords <- formatted_data[, c("Cell.X.Position", 
                                       "Cell.Y.Position")]
      db <- dbscan::dbscan(cell_cords, eps = radius, minPts = min_neighborhood_size)
      formatted_data$Cluster <- factor(db$cluster + 1)
    }
    else if (method == "rphenograph") {
      stop("This option is not available for this version yet! Check dev version for this function!")
    }
    else {
      stop("Please select a valid clustering method, current options: dbscan")
    }
  }
  cells_in_clusters <- formatted_data[stats::complete.cases(formatted_data), 
  ]
  cells_not_in_clusters <- formatted_data[!stats::complete.cases(formatted_data), 
  ]
  number_of_clusters <- length(unique(cells_in_clusters$Cluster))
  if (plot_clusters) {
    label_location <- vector()
    for (Clusternumber in seq_len(number_of_clusters)) {
      cells_in_Cluster <- cells_in_clusters[cells_in_clusters$Cluster == 
                                              Clusternumber, ]
      minX <- min(cells_in_Cluster$Cell.X.Position, na.rm = TRUE)
      maxX <- max(cells_in_Cluster$Cell.X.Position, na.rm = TRUE)
      minY <- min(cells_in_Cluster$Cell.Y.Position, na.rm = TRUE)
      maxY <- max(cells_in_Cluster$Cell.Y.Position, na.rm = TRUE)
      averageX <- (minX + maxX)/2
      averageY <- (minY + maxY)/2
      label_location <- rbind(label_location, c(Clusternumber, 
                                                averageX, averageY))
    }
    label_location <- as.data.frame(label_location)
    colnames(label_location) <- c("Cluster", "Xpos", "Ypos")
    cluster_colours <- dittoSeq::dittoColors()[seq_len(number_of_clusters)]
    q <- ggplot(cells_in_clusters, aes(x = Cell.X.Position, 
                                       y = Cell.Y.Position))
    q <- q + geom_point(aes(color = Cluster))
    q <- q + geom_text(data = label_location, aes(x = Xpos, 
                                                  y = Ypos, label = Cluster))
    q <- q + scale_color_manual(values = cluster_colours)
    if (dim(cells_not_in_clusters)[1] != 0) {
      q <- q + geom_point(data = cells_not_in_clusters, 
                          colour = "black")
    }
    q <- q + xlab("Cell.X.Position") + ylab("Cell.Y.Position") + 
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), legend.position = "none")
    methods::show(q)
  }
  formatted_data_with_clusters <- formatted_data
  formatted_data_with_clusters$Cluster <- paste0("Cluster_", 
                                                 as.character(formatted_data_with_clusters$Cluster))
  formatted_data_with_clusters$Cluster[formatted_data_with_clusters$Cluster == 
                                         "Cluster_NA"] <- "Free_cell"
  formatted_data_with_clusters$Cluster[formatted_data_with_clusters$Cluster == 
                                         "Cluster_Free_cell"] <- "Free_cell"
  colnames(formatted_data_with_clusters)[which(colnames(formatted_data_with_clusters) == 
                                                 "Cluster")] <- "Neighborhood"
  if (is.null(SummarizedExperiment::colData(spe_object)$Cell.ID)) {
    SummarizedExperiment::colData(spe_object)$Cell.ID <- rownames(SummarizedExperiment::colData(spe_object))
  }
  SummarizedExperiment::colData(spe_object) <- methods::as(merge(data.frame(SummarizedExperiment::colData(spe_object)), 
                                                                 formatted_data_with_clusters[, c("Cell.ID", "Neighborhood")], 
                                                                 by = "Cell.ID", all.x = TRUE), "DFrame")
  return(spe_object)
}

for (target in c("p53", "CD3", "CD8", "CD20")) {
  spe <- F4_B2
  reference <- "p53"
  cat(target, ": ", 
  (mean(na.omit(
  cell_density_in_neighborhood(spe, reference_cell =reference, target_cell = target, radius = 200)
  )))/(unlist(table(spe$Phenotype))[target]/(pi*1500^2)),
  "\n")
  }

test1 <- cell_density_in_neighborhood(F4_B2, "CD3", "CD3", 100)
test2 <- cell_density_in_neighborhood(F4_B2, "CD3", "CD8", 100)
test3 <- cell_density_in_neighborhood(F4_B2, "CD3", "CD20", 100)
test4 <- cell_density_in_neighborhood(F4_B2, "CD3", "p53", 100)

