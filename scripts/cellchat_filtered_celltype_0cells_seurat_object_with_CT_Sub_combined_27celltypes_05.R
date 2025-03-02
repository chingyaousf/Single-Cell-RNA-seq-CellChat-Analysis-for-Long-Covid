# Load the required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)


## comparing IFN-II pathway in Longcovid and control
# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("IFN-II")

# Check which datasets have significant communication for the pathway
datasets_to_plot <- list()
for (i in 1:length(object.list)) {
  if (pathways.show %in% object.list[[i]]@netP$pathways) {
    datasets_to_plot[[names(object.list)[i]]] <- object.list[[i]]
  }
}

# If no datasets have significant communication, print a message
if (length(datasets_to_plot) == 0) {
  stop("No significant communication found for the IFN-II pathway in any dataset.")
}

# Generate Circle plots only for datasets with significant communication
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/IFN-II_circle_comparison_updated.png",
    width = 8 * length(datasets_to_plot), height = 8, units = "in", res = 300)

# Set up the plotting layout
par(mfrow = c(1, length(datasets_to_plot)), xpd = TRUE)
for (i in names(datasets_to_plot)) {
  netVisual_aggregate(datasets_to_plot[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, i))
}
dev.off()


## IFN-II pathway only in Longcovid
# Load the saved CellChat object for LongCovid
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")

# Check the available pathways
#pathways_available <- longcovid_cellchat@netP$pathways
#print(pathways_available)

# Define the signaling pathway to visualize (IFN-II)
pathways.show <- c("IFN-II")

# Verify if IFN-II is a significant pathway in the LongCovid dataset
#if (!pathways.show %in% pathways_available) {
#  stop("IFN-II pathway is not significant in the LongCovid dataset.")
#}

# Aggregate the communication network
longcovid_cellchat <- aggregateNet(longcovid_cellchat)

# Optional: Hierarchy plot (uncomment if needed)
# Define the receiver nodes (e.g., immune cells or specific clusters)
# vertex.receiver <- seq(1, 4) # Adjust indices as needed
# netVisual_aggregate(longcovid_cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Circle plot for the IFN-II pathway
par(mfrow = c(1, 1))  # Single plot layout
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/IFN-II_circle_LongCovid.png",
    width = 8, height = 8, units = "in", res = 300)
netVisual_aggregate(longcovid_cellchat, signaling = pathways.show, layout = "circle")
dev.off()




## ICAM_circle_comparison, (not) showing interactive windows before save plot
# Load the required packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# Generate Circle plots for ICAM pathway and save them as a combined plot
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_circle_comparison_updated.png",
    width = 16, height = 8, units = "in", res = 300)

# Set up the plotting layout and create Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()


## ICAM_circle_comparison, showing interactive windows before save plot
# Load the required packages
library(Seurat)
library(CellChat)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# ---- STEP 1: Review the plots interactively ----
# Set up the plotting layout and display Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# ---- STEP 2: Save the plots after review ----
# Save the Circle plots to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_circle_comparison_updated.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Circle plots for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()




## ICAM_circle_comparison between Longcovid and control, for (CD4.Tfh) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD4.Tfh) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions between the CD4.Tfh cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD4.Tfh as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway (can be modified)

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD4.Tfh"] <- "CD4.Tfh"  # Highlight CD4.Tfh cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.Tfh_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()




## ICAM_circle_comparison between Longcovid and control, for (CD8.TE) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD8.TE) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions between the CD8.TE cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD8.TE as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway (can be modified)

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD8.TE"] <- "CD8.TE"  # Highlight CD8.TE cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD8.TE_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()




## CypA_circle_comparison between Longcovid and control, for (CD4.Tfh) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD8.TE) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (CypA) signaling pathway 
# to visualize the interactions between the CD8.TE cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD8.TE as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("CypA")  # Updated to the CypA pathway

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD8.TE"] <- "CD8.TE"  # Highlight CD8.TE cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/CypA_chord_CD8.TE_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()



## MHC-II_circle_comparison between Longcovid and control, for (CD8.TE) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD8.TE) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (MHC-II) signaling pathway 
# to visualize the interactions between the CD8.TE cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD8.TE as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("MHC-II")  # Updated to the MHC-II pathway

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD8.TE"] <- "CD8.TE"  # Highlight CD8.TE cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/MHC-II_chord_CD8.TE_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()



## MHC-II_circle_comparison between Longcovid and control, for (CD4.Tfh) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD4.Tfh) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (MHC-II) signaling pathway 
# to visualize the interactions between the CD4.Tfh cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD4.Tfh as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("MHC-II")  # Updated to the MHC-II pathway

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD4.Tfh"] <- "CD4.Tfh"  # Highlight CD4.Tfh cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/MHC-II_chord_CD4.Tfh_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()



## MHC-II_circle_comparison between Longcovid and control, for (CD4.Tfh) Cell Type (Receiving Signals) Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD4.Tfh) Cell Type (Receiving Signals) from All Other Cell Types
# Description: This script generates chord diagrams for the (MHC-II) signaling pathway 
# to visualize the signals received by the CD4.Tfh cell type from all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD4.Tfh as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("MHC-II")  # Updated to the MHC-II pathway

# Define the cell type of interest (receiver) and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD4.Tfh"] <- "CD4.Tfh"  # Highlight CD4.Tfh cell type as receiver

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       sources.use = NULL,        # Do not filter by sources
                       targets.use = "CD4.Tfh",   # Focus on signals received by CD4.Tfh
                       title.name = paste0(pathways.show, " signaling received by CD4.Tfh - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/MHC-II_chord_received_CD4.Tfh_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       sources.use = NULL,        # Do not filter by sources
                       targets.use = "CD4.Tfh",   # Focus on signals received by CD4.Tfh
                       title.name = paste0(pathways.show, " signaling received by CD4.Tfh - ", names(object.list)[i]))
}

dev.off()






## ICAM_circle_CD8.TE_target_Bcells_source
# Title: Generate Circle Plots for ICAM Signaling Between CD8.TE (target) and Selected B Cell Types (sources)
# Description: This script generates circle plots for the (ICAM) signaling pathway, 
# focusing only on interactions where CD8.TE is the target and the selected B cell types 
# ("B_naive", "B_exhausted", "B_immature", "B_malignant", "B_memory") are the sources.

# Load the required packages
library(Seurat)
library(CellChat)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")

# Define the sources (B cell types) and targets (CD4.Tfh)
sources.use <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", "B_memory")
targets.use <- c("CD8.TE")

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# ---- STEP 1: Review the plots interactively ----
# Set up the plotting layout and display Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      sources.use = sources.use,    # Focus only on B cell types as sources
                      targets.use = targets.use,    # Focus only on CD4.Tfh as the target
                      remove.isolate = FALSE)
}

# ---- STEP 2: Save the plots after review ----
# Save the Circle plots to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_circle_CD8.TE_target_Bcells_source.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Circle plots for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      sources.use = sources.use,    # Focus only on B cell types as sources
                      targets.use = targets.use,    # Focus only on CD8.TE as the target
                      remove.isolate = FALSE)
}

dev.off()



## MHC-II_circle_CD4.Tfh_target_Bcells_source
# Title: Generate Circle Plots for MHC-II Signaling Between CD4.Tfh (target) and Selected B Cell Types (sources)
# Description: This script generates circle plots for the (MHC-II) signaling pathway, 
# focusing only on interactions where CD4.Tfh is the target and the selected B cell types 
# ("B_naive", "B_exhausted", "B_immature", "B_malignant", "B_memory") are the sources.

# Load the required packages
library(Seurat)
library(CellChat)
library(ggplot2)

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("MHC-II")  # Updated to MHC-II signaling pathway

# Define the sources (B cell types) and targets (CD4.Tfh)
sources.use <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", "B_memory")
targets.use <- c("CD4.Tfh")  # Updated target to CD4.Tfh

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# ---- STEP 1: Review the plots interactively ----
# Set up the plotting layout and display Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      sources.use = sources.use,    # Focus only on B cell types as sources
                      targets.use = targets.use,    # Focus only on CD4.Tfh as the target
                      remove.isolate = FALSE)
}

# ---- STEP 2: Save the plots after review ----
# Save the Circle plots to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/MHC-II_circle_CD4.Tfh_target_Bcells_source.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Circle plots for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      sources.use = sources.use,    # Focus only on B cell types as sources
                      targets.use = targets.use,    # Focus only on CD4.Tfh as the target
                      remove.isolate = FALSE)
}

dev.off()






## ICAM_circle_comparison between Longcovid and control, for (CD4.IL22) Cell Type
# Title: Generate Chord Diagrams for (CD4.IL22) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions of the CD4.IL22 cell type with all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD4.IL22 as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)


# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway

# Define the cell type of interest: CD4.IL22
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD4.IL22"] <- "CD4.IL22"  # Highlight CD4.IL22 cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD4.IL22- ", names(object.list)[i]))
}

# ---- STEP 2: Save the Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.IL22_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD4.IL22- ", names(object.list)[i]))
}

dev.off()




## ICAM_circle_comparison between Longcovid and control, for (gdT) Cell Type
# Title: Generate Chord Diagrams for (gdT) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions of the gdT cell type with all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight gdT as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway

# Define the cell type of interest: gdT
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["gdT"] <- "gdT"  # Highlight gdT cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network gdT- ", names(object.list)[i]))
}

# ---- STEP 2: Save the Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_gdT_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network gdT- ", names(object.list)[i]))
}

dev.off()





## TGFb_circle_comparison, showing interactive windows before saving the plot
# Load the required packages
library(Seurat)
library(CellChat)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("TGFb")  # Update to TGFb signaling pathway

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# ---- STEP 1: Review the plots interactively ----
# Set up the plotting layout and display Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# ---- STEP 2: Save the plots after review ----
# Save the Circle plots to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/TGFb_circle_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Circle plots for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()





## TGFb_circle_comparison between Longcovid and control, for (CD14_mono) Cell Type Interacting with All Other Cell Types
# Title: Generate Chord Diagrams for (CD14_mono) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (TGFb) signaling pathway 
# to visualize the interactions between the CD14_mono cell type and all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD14_mono as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("TGFb")  # Update to TGFb signaling pathway

# Define the cell type of interest and group all others
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD14_mono"] <- "CD14_mono"  # Highlight CD14_mono cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/TGFb_chord_CD14_mono_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()





## JAM_circle_comparison, showing interactive windows before saving the plot
## control no JAM signal, only in LongCovid, so this code error message, no JAM signal in control
# Load the required packages
library(Seurat)
library(CellChat)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("JAM")  # Update to JAM signaling pathway

# Get the maximum weight for consistent edge weight scaling across datasets
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

# ---- STEP 1: Review the plots interactively ----
# Set up the plotting layout and display Circle plots for each dataset
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# ---- STEP 2: Save the plots after review ----
# Save the Circle plots to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/JAM_circle_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Circle plots for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for side-by-side plots
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()




## JAM pathway only in Longcovid
# Load the saved CellChat object for LongCovid
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")

# Check the available pathways
#pathways_available <- longcovid_cellchat@netP$pathways
#print(pathways_available)

# Define the signaling pathway to visualize (IFN-II)
pathways.show <- c("JAM")

# Verify if IFN-II is a significant pathway in the LongCovid dataset
#if (!pathways.show %in% pathways_available) {
#  stop("IFN-II pathway is not significant in the LongCovid dataset.")
#}

# Aggregate the communication network
longcovid_cellchat <- aggregateNet(longcovid_cellchat)

# Optional: Hierarchy plot (uncomment if needed)
# Define the receiver nodes (e.g., immune cells or specific clusters)
# vertex.receiver <- seq(1, 4) # Adjust indices as needed
# netVisual_aggregate(longcovid_cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Circle plot for the IFN-II pathway
par(mfrow = c(1, 1))  # Single plot layout
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/JAM_circle_LongCovid.png",
    width = 8, height = 8, units = "in", res = 300)
netVisual_aggregate(longcovid_cellchat, signaling = pathways.show, layout = "circle")
dev.off()






## ICAM_circle_comparison between Longcovid and control, for (CD4.Naive) Cell Type
# Title: Generate Chord Diagrams for (CD4.Naive) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions of the CD4.Naive cell type with all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD4.Naive as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway

# Define the cell type of interest: CD4.Naive
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD4.Naive"] <- "CD4.Naive"  # Highlight CD4.Naive cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD4.Naive - ", names(object.list)[i]))
}

# ---- STEP 2: Save the Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.Naive_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD4.Naive - ", names(object.list)[i]))
}

dev.off()




## ICAM_circle_comparison between Longcovid and control, for (CD8.EM) Cell Type
# Title: Generate Chord Diagrams for (CD8.EM) Cell Type Interacting with All Other Cell Types
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions of the CD8.EM cell type with all other cell types 
# in both Control and LongCovid datasets. The diagrams highlight CD8.EM as a distinct group.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway

# Define the cell type of interest: CD8.EM
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)
group.cellType["CD8.EM"] <- "CD8.EM"  # Highlight CD8.EM cell type

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD8.EM - ", names(object.list)[i]))
}

# ---- STEP 2: Save the Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD8.EM_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network CD8.EM - ", names(object.list)[i]))
}

dev.off()





# Title: Generate Chord Diagrams for (ICAM) Signaling in CD4.Tfh, CD4.IL22, CD8.TE, and gdT
# Description: This script generates chord diagrams for the (ICAM) signaling pathway 
# to visualize the interactions between specific cell types (CD4.Tfh, CD4.IL22, CD8.TE, and gdT) 
# and all other cell types in both Control and LongCovid datasets.

# Load the required packages
library(Seurat)
library(CellChat)
library(circlize)
library(ggplot2)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_27celltypes.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_27celltypes.rds")


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Create a list of CellChat objects for comparison
object.list <- list(Control = control_cellchat, LongCovid = longcovid_cellchat)

# Define the signaling pathway to visualize
pathways.show <- c("ICAM")  # Example pathway

# Define the cell types of interest and group all others as "Other"
group.cellType <- rep("Other", length(levels(object.list[[1]]@idents)))
names(group.cellType) <- levels(object.list[[1]]@idents)

# Highlight the selected cell types
highlighted_cells <- c("CD4.Tfh", "CD4.IL22", "CD8.TE", "gdT")
group.cellType[highlighted_cells] <- highlighted_cells  

# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.Tfh_CD4.IL22_CD8.TE_gdT_comparison.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

dev.off()


## increase label text
# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]),
                       lab.cex = 1.0     # Increase label size
  )
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.Tfh_CD4.IL22_CD8.TE_gdT_comparison_02.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]),
                       lab.cex = 1.0     # Increase label size
                       
  )
}

dev.off()


## increase title text, but overlap with small one
# ---- STEP 1: Display Chord Diagrams interactively ----
# Set up the plotting layout for interactive viewing
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid

for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = NULL,  # Remove the built-in title
                       lab.cex = 1.0  # Increase label size
  )
  title(main = paste0(pathways.show, " signaling network - ", names(object.list)[i]), 
        cex.main = 2.5, font.main = 2)  # Add title separately with size and bold
}

# ---- STEP 2: Save Chord Diagrams after review ----
# Save the Chord Diagrams to a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/redo/ICAM_chord_CD4.Tfh_CD4.IL22_CD8.TE_gdT_comparison_03.png",
    width = 16, height = 8, units = "in", res = 300)

# Re-create the Chord Diagrams for saving
par(mfrow = c(1, 2), xpd = TRUE)  # Side-by-side plots for Control and LongCovid

for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], 
                       signaling = pathways.show, 
                       group = group.cellType, 
                       title.name = NULL,  # Remove the built-in title
                       lab.cex = 1.0  # Increase label size
  )
  title(main = paste0(pathways.show, " signaling network - ", names(object.list)[i]), 
        cex.main = 2.5, font.main = 2)  # Add title separately with size and bold
}

dev.off()

