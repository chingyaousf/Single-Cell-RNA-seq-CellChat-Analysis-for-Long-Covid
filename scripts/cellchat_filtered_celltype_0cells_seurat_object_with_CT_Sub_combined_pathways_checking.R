# Load required packages
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)

# Load the CellChat object
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")

# Get all signaling pathways identified in the CellChat object
signaling_pathways <- longcovid_cellchat@netP$pathways

# Display the pathways
print(signaling_pathways)


# Load the saved CellChat object
seurat.cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/seurat_cellchat_combined.RDS")

# Get all signaling pathways identified in the CellChat object
signaling_pathways <- seurat.cellchat@netP$pathways

# Display the pathways
print(signaling_pathways)
