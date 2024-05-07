#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexHeatmap)


mIHC_counts1 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_1_dataframe.csv")
mIHC_counts2 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_2_dataframe.csv")
mIHC_counts3 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_3_dataframe.csv")
mIHC_counts4 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_4_dataframe.csv")
mIHC_counts5 <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/cell_count_5_dataframe.csv")


mIHC_counts <- rbind(mIHC_counts1, mIHC_counts2, mIHC_counts3, mIHC_counts4, mIHC_counts5)

#prel_mIHC_counts <- read.csv("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/cell_counts/prel_cell_count_dataframe.csv")
#rownames(prel_mIHC_counts) = prel_mIHC_counts$CORE_ID
rownames(mIHC_counts) = mIHC_counts$CORE_ID



sIHC_counts <- read.csv2("/home/isc/Spatial_immune_env/data_from_suze/data/supplData_withimages.csv")

# Reading cell counts from TMArQ
tmarq_counts <- read.table("/home/isc/Spatial_immune_env/vectra/segmentation_and_phenotyping/counts_from_TMArQ/segmentcell_multiplex.txt", header = TRUE)

tmarq_counts$uid <- sapply(strsplit(tmarq_counts$file, '[][]'), function(element) {if (!is.na(element[2])) 
                                                                                    {paste0("BLOCK_", 
                                                                                      strsplit(element[2], ",")[[1]][1], 
                                                                                      "|",
                                                                                      strsplit(element[2], ",")[[1]][3],
                                                                                      strsplit(element[2], ",")[[1]][2])}
                                                                                   else {NA}})

tmarq_counts$marker <- sapply(strsplit(tmarq_counts$file, '_'), function(element) element[5])

for (element in rownames(mIHC_counts)) {
  
cat ("\n", element, ": ", abs(prel_mIHC_counts[element, "PAN.CK"] - mIHC_counts[element, "PAN.CK"]))  
  
}

#
# PRELIMINARY VS FINAL PHENOTYPE
#

# plot(x=mIHC_counts$PAN.CK, 
#      y=prel_mIHC_counts$PAN.CK,
#      xlab ="Final counts counts",
#      ylab ="Preliminary counts",
#      main = "PAN-CK COUNTS",
#      pch=16)
# abline(1,1)
# 
# plot(x=mIHC_counts$CD8, 
#      y=prel_mIHC_counts$CD8,
#      xlab ="Final counts counts",
#      ylab ="Preliminary counts",
#      main = "CD8 COUNTS",
#      pch=16)
# abline(1,1)
# 
# plot(x=mIHC_counts$CD4, 
#      y=prel_mIHC_counts$CD4,
#      xlab ="Final counts counts",
#      ylab ="Preliminary counts",
#      main = "CD4 COUNTS",
#      pch=16)
# abline(1,1)
# 
# plot(x=mIHC_counts$CD20, 
#      y=prel_mIHC_counts$CD20,
#      xlab ="Final counts counts",
#      ylab ="Preliminary counts",
#      main = "CD20 COUNTS",
#      pch=16)
# abline(1,1)

#
# TOTAL COUNTS VS TMArQ ESTIMATED TOTAL COUNT
#

tot_counts <- apply(mIHC_counts[,c(seq(6,13))], FUN=sum, MARGIN=1)
filt1 <- tmarq_counts[tmarq_counts$marker == "DAPI", ]
rownames(filt1) <- filt1$uid

data <- data.frame(
  Vectra_Polaris_counts = tot_counts,
  TMArQ_counts = filt1[names(tot_counts), "hema_cells"]
)


# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts",) +
  title("All") +  xlim(0,15000) +
  ylim(0,15000) +
  theme_classic()


#
# FINAL PHENOTYPE VS TMArQ ESTIMATED FROM SAME SLICE
#

# Getting info of interest
filt1 <- tmarq_counts[tmarq_counts$marker == "CD8", ]
rownames(filt1) <- filt1$uid

data = data.frame(Vectra_Polaris_counts=mIHC_counts$CD8 + mIHC_counts$CD8_FOXP3,
                  TMArQ_counts=filt1[mIHC_counts$CORE_ID, "dab_cells10"])

# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts") +
  xlim(0,4500) +
  ylim(0,4500) +
  theme_classic()


# Getting info of interest
filt1 <- tmarq_counts[tmarq_counts$marker == "CD4", ]
rownames(filt1) <- filt1$uid

data = data.frame(Vectra_Polaris_counts=mIHC_counts$CD4 + mIHC_counts$CD4_FOXP3,
                  TMArQ_counts=filt1[mIHC_counts$CORE_ID, "dab_cells10"])

# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts") +
  xlim(0,4500) +
  ylim(0,4500) +
  theme_classic()


filt1 <- tmarq_counts[tmarq_counts$marker == "CD20", ]
rownames(filt1) <- filt1$uid


data = data.frame(Vectra_Polaris_counts=mIHC_counts$CD20,
                  TMArQ_counts=filt1[mIHC_counts$CORE_ID, "dab_cells10"])

# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts") +
  xlim(0,4500) +
  ylim(0,4500) +
  theme_classic()

filt1 <- tmarq_counts[tmarq_counts$marker == "pan-CK", ]
rownames(filt1) <- filt1$uid


data = data.frame(Vectra_Polaris_counts=mIHC_counts$PAN.CK,
                  TMArQ_counts=filt1[mIHC_counts$CORE_ID, "dab_cells10"])

# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts") +
  xlim(0,8000) +
  ylim(0,8000) +
  theme_classic()


filt1 <- tmarq_counts[tmarq_counts$marker == "FOXP3", ]
rownames(filt1) <- filt1$uid


data = data.frame(Vectra_Polaris_counts=mIHC_counts$CD8_FOXP3 + mIHC_counts$CD4_FOXP3,
                  TMArQ_counts=filt1[mIHC_counts$CORE_ID, "dab_cells10"])

# Create the ggplot
ggplot(data, aes(x = Vectra_Polaris_counts, y = TMArQ_counts)) +
  geom_point(shape = 16) +  # scatter plot
  geom_abline(slope = 1, intercept = 0, col="darkgreen") +  # line y = x
  labs(
    x = "Phenoimager counts",
    y = "TMArQ counts") +
  xlim(0,1000) +
  ylim(0,1000) +
  theme_classic()


# #
# # COMPARING RESULTS FROM sIHC from different slices
# #
# 
# 
# #Usiing uid as rowname for sihc counts
# rownames(sIHC_counts) <- sIHC_counts$uid
# 
# # Plotting CD8 counts and log counts
# plot(x=log(mIHC_counts$CD8), 
#      y=log(sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD8"]),
#      xlab ="Vectra Polaris log counts",
#      ylab ="TMArQ log counts",
#      main = "log CD8 counts",
#      pch=16)
# 
# abline(a=0, b=1, lwd=6, col="green")
# 
# # Plotting CD8 counts and log counts
# plot(x=mIHC_counts$CD8, 
#      y=sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD8"],
#      xlab ="Vectra Polaris counts",
#      ylab ="TMArQ counts",
#      main = "CD8 counts",
#      pch=16)
# abline(a=0, b=1, lwd=6, col="green")
# 
# 
# 
# # Plotting CD4 counts and log counts
# plot(x=log(mIHC_counts$CD4), 
#      y=log(sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD4"]),
#      xlab ="Vectra Polaris log counts",
#      ylab ="TMArQ log counts",
#      main = "log CD4 counts",
#      pch=16)
# abline(a=0, b=1, lwd=6, col="green")
# 
# # Plotting CD4 counts and log counts
# plot(x=mIHC_counts$CD4, 
#      y=sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD4"],
#      xlab ="Vectra Polaris counts",
#      ylab ="TMArQ counts",
#      main = "CD4 counts",
#      pch=16)
# abline(a=0, b=1, lwd=6, col="green")
# 
# 
# 
# # Plotting CD20 counts and log counts
# plot(x=log(mIHC_counts$CD20), 
#      y=log(sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD20"]),
#      xlab ="Vectra Polaris log counts",
#      ylab ="TMArQ log counts",
#      main = "log CD20 counts",
#      pch=16)
# abline(a=0, b=1, lwd=6, col="green")
# 
# # Plotting CD4 counts and log counts
# plot(x=mIHC_counts$CD20, 
#      y=sIHC_counts[mIHC_counts$TMArQ_CORE_ID, "CD20"],
#      xlab ="Vectra Polaris counts",
#      ylab ="TMArQ counts",
#      main = "CD20 counts",
#      pch=16)
# 
# abline(a=0, b=1, lwd=6, col="green")



