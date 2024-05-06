#!/usr/bin/env Rscript

# Install and load the required packages
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  install.packages("ComplexHeatmap")

library(ComplexHeatmap)

if (!requireNamespace("Rtsne", quietly = TRUE))
  install.packages("Rtsne")

library(Rtsne)

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(ggplot2)

# Load the intensity data of an example phenotyped mIHC core
data <- read.csv("/media/isc/Data1/Processed_cores/SCAN-B_TNBC_TMA_1A/20220929_1A_Core[1,1,D]_[10938,41418]_component_data.tif_intensity.csv", header = TRUE)

# Determining and plotting tSNE
my_tsne <- Rtsne(data[,seq(5,10)])


# Set specific colors for each phenotype
unique_phenotypes <- sort(unique(data$Phenotype))
num_phenotypes <- length(unique_phenotypes)
phenotype_colors <- c("yellow", "green", "darkgreen", "purple", "red", "darkred", "grey", "pink")


#Plotting scatterplot of the core

ggplot(data=data, aes(x=X_coor, y=Y_coor)) + 
  geom_point(aes(col=Phenotype)) +
  scale_color_manual(values = setNames(phenotype_colors, unique_phenotypes)) +  # Set specific colors manually
  theme_classic()

# Plot tSNE
ggplot(data = as.data.frame(my_tsne$Y), aes(x = V1, y = V2)) +
  geom_point(aes(col = data$Phenotype)) +
  scale_color_manual(values = setNames(phenotype_colors, unique_phenotypes)) +  # Set specific colors manually
  theme_classic() +
  labs(col = "Determined phenotype", x = "tSNE 1", y = "tSNE 2")


# Plot heatmap
phenotype_annotation <- rowAnnotation(Phenotype = as.factor(data$Phenotype),
                                      col = list(Phenotype=setNames(phenotype_colors, unique_phenotypes)))

order_of_rows <- order(data$Phenotype)

Heatmap(log(data[,seq(5,10)] + 1e-50),
        name = "log Intensity",
        row_order = order_of_rows,
        cluster_columns = FALSE,
        left_annotation = phenotype_annotation)
