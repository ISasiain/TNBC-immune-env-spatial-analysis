rm(list = ls())
setwd("/home/isc/Spatial_immune_env")


# Installing packages if requitred
if (!require("survival", quietly = TRUE))
  install.packages("survivalr")

if (!require("survminer", quietly = TRUE))
  install.packages("survminer")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("SPIAT", quietly = TRUE))
  install.packages("SPIAT")


library(ggplot2)
library(survival)
library(survminer)
library(SPIAT)

#Reading annotation file and generating survival plots
annotation <- read.csv2("data_from_suze/data/supplData_withimages.csv")
rownames(annotation) <- annotation$uid



#
# Read multiple marker statistics file
#

multiple_marker_stats <- readRDS("CD3_CD8_spatial_stats.RData")
ref_cell_type <- c("CD3")
target_cell_type <- c("CD20")




#
# Plot the effect of AMD in survival
#

AMD_Median <- unlist(lapply(multiple_marker_stats, function(item) item$AMD[1, "Median"]))
names_AMD <- names(unlist(lapply(multiple_marker_stats, function(item) item$AMD[1, "Median"])))

plot(sort(AMD_Median))

analyse_AMD <- as.data.frame(cbind(
  
  "name"=names_AMD,
  "time"=as.numeric(annotation[names_AMD, "time"]),
  "event"= annotation[names_AMD, "event"],
  "Median_AMD"= cut(AMD_Median,
                  breaks = c(0, 20, 30, 40, 100, 200),
                  labels = c("0-20", "20-30","30-40" , "40-100", "100-200"),
                  include_lowest = TRUE)
  
))

analyse_AMD$time <- as.numeric(analyse_AMD$time)
analyse_AMD$event <- as.numeric(analyse_AMD$event)
analyse_AMD$Median_AMD <- as.factor(analyse_AMD$Median_AMD)


fit <- survfit(Surv(time, event) ~ Median_AMD, data = analyse_AMD)

ggsurvplot(fit, 
           data=analyse_AMD,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)")



#
# Plot the effect of CIN in survival
#

CIN <- unlist(lapply(multiple_marker_stats, function(item) item$CIN))
names_CIN <- names(unlist(lapply(multiple_marker_stats, function(item) item$CIN)))

hist(sort(CIN))
sort(CIN)

analyse_CIN <- as.data.frame(cbind(
  
  "name"=names_CIN,
  "time"=as.numeric(annotation[names_CIN, "time"]),
  "event"= annotation[names_CIN, "event"],
  "CIN"= cut(CIN,
                 breaks = c(0, 3, 5, 10, 20, 30, 100),
                 labels = c("0-3", "3-5", "5-10", "10-20", "20-30", "30-100"),
                 include_lowest = TRUE)
  
))

analyse_CIN$time <- as.numeric(analyse_CIN$time)
analyse_CIN$event <- as.numeric(analyse_CIN$event)
analyse_CIN$CIN <- as.factor(analyse_CIN$CIN)


fit <- survfit(Surv(time, event) ~ CIN, data = analyse_CIN)

ggsurvplot(fit, 
           data=analyse_CIN,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)", 
           label.title="CIN (%)",
           legend.labs=c("0-3", "3-5", "5-10", "10-20", "20-30", "30-100")
)




#
# Plot the effect of attraction in survival
#

obs_vs_expected <- sapply(multiple_marker_stats, function(item) unname(item$"CIN" / item$"Cell_Proportions"[target_cell_type]))

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]


attraction <- unlist(obs_vs_expected)
names_attraction <- names(attraction)

plot(sort(attraction), ylim=c(0,2))

analyse_attraction<- as.data.frame(cbind(
  
  "name"=names_attraction,
  "time"=as.numeric(annotation[names_attraction, "time"]),
  "event"= annotation[names_attraction, "event"],
  "attraction"= cut(attraction,
                    breaks = c(0, 0.5, 1, 1.5, 2, 3, 5),
                    labels = c("0-0.5", "0.5-1", "1-1.5", "1.5-2", "2-3", "3-5"),
                    include_lowest = TRUE)
  
  
))



analyse_attraction$time <- as.numeric(analyse_attraction$time)
analyse_attraction$event <- as.numeric(analyse_attraction$event)


fit <- survfit(Surv(time, event) ~ attraction, data = analyse_attraction)


ggsurvplot(fit, 
           data=analyse_attraction,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)")


#
# Plot the effect of clustering
#

K_curve <- sapply(multiple_marker_stats, function(item) item$"AUC")
AUC <- sapply(K_curve, function(item) AUC_of_cross_function(item))


names_AUC <- names(AUC)


analyse_AUC <- as.data.frame(cbind(
  
  "name"=names_AUC,
  "time"=as.numeric(annotation[names_AUC, "time"]),
  "event"= annotation[names_AUC, "event"],
  "AUC"= cut(AUC,
                    breaks = c(-100, -0.1, 0.1, 100),
                    labels = c("Ind. clusters", "Random", "Clustered together"),
                    include_lowest = TRUE)
  
  
))

analyse_AUC$time <- as.numeric(analyse_AUC$time)
analyse_AUC$event <- as.numeric(analyse_AUC$event)

fit <- survfit(Surv(time, event) ~ AUC, data = analyse_AUC)


ggsurvplot(fit, 
           data=analyse_AUC,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)")


