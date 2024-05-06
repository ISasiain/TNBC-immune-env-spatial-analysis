#!/usr/bin/env Rscript

# Install and load the required packages

if (!requireNamespace("ComplexHeatamap", quietly = TRUE))
    install.packages("ComplexHeatmap")

library(ComplexHeatmap)

if (!requireNamespace("Rtsne", quietly = TRUE))
    install.packages("Rtsne")

library(Rtsne)


# Load the intensity data of an example phenotyped mIHC core

data <- read.csv("/media/isc/Data1/Processed_cores/SCAN-B_TNBC_TMA_1A/20220929_1A_Core[1,1,D]_[10938,41418]_component_data.tif_intensity.csv")
