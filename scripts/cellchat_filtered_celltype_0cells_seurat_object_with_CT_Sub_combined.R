# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(CellChat)
#library(ggalluvial)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined.rds")

# Inspect the Seurat object
head(seurat_object@meta.data)
seurat_object
unique(seurat_object@meta.data$CT_Sub_combined)

# Set up CellChat database
CellChatDB <- CellChatDB.human

# Set cell identities to the "CT_Sub_combined" column and specify the default assay
Idents(seurat_object) <- seurat_object$CT_Sub_combined
DefaultAssay(seurat_object) <- "RNA"

# Create CellChat object
seurat.cellchat <- createCellChat(object = seurat_object, group.by = "CT_Sub_combined")

# Drop unused factor levels from cell identities
seurat.cellchat@idents <- droplevels(seurat.cellchat@idents)

# Check if the data matrix matches between Seurat and CellChat objects
summary(seurat.cellchat@data@i == seurat_object@assays[["RNA"]]@data@i)['FALSE'] # NA if no false

# Set the CellChatDB for signaling
seurat.cellchat@DB <- CellChatDB

# Pre-processing
seurat.cellchat <- subsetData(seurat.cellchat) # Subset data based on the CellChatDB
seurat.cellchat <- identifyOverExpressedGenes(seurat.cellchat) # Identify overexpressed genes
seurat.cellchat <- identifyOverExpressedInteractions(seurat.cellchat) # Identify overexpressed interactions

# Compute the communication probability
seurat.cellchat <- computeCommunProb(seurat.cellchat)

# Infer cellular communication network
seurat.cellchat <- filterCommunication(seurat.cellchat, min.cells = 10)
communication.df <- subsetCommunication(seurat.cellchat)

# Inspect the Data
print(head(communication.df)) # Display the first few rows of the data frame to inspect

# Save the Data Frame
write.csv(communication.df, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/communication_interactions_combined.csv", row.names = FALSE)

# Infer the cell-cell communication at a signaling pathway level
seurat.cellchat <- computeCommunProbPathway(seurat.cellchat)
saveRDS(seurat.cellchat, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/seurat_cellchat_combined.RDS")

# Load the saved CellChat object
seurat.cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/seurat_cellchat_combined.RDS")

head(seurat.cellchat@idents)
slotNames(seurat.cellchat)
str(seurat.cellchat)
unique(seurat.cellchat@meta$disease)


# Calculate the aggregated cell-cell communication network
seurat.cellchat <- aggregateNet(seurat.cellchat)

# Visualize the number of interactions with weights/strengths
matrix <- as.data.frame(seurat.cellchat@net$count)
matrix[matrix == 0] <- 0.000000001
matrix <- data.matrix(matrix)

# Visualization of interaction numbers
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_number_combined.png", 
    width = 10,
    height = 10, 
    units = "in",
    res = 500)
netVisual_circle(seurat.cellchat@net$count,
                 top = 0.1,
                 arrow.width = 2,
                 arrow.size = 0.5,
                 edge.width.max = 10, 
                 weight.scale = TRUE, 
                 label.edge = FALSE)
dev.off()

# Visualization of interaction strengths
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_strength_combined.png", 
    width = 10,
    height = 10, 
    units = "in",
    res = 500)
netVisual_circle(seurat.cellchat@net$weight,
                 top = 0.1,
                 arrow.width = 2,
                 arrow.size = 0.5,
                 edge.width.max = 10, 
                 weight.scale = TRUE, 
                 label.edge = FALSE)
dev.off()

# Save plots of all inferred networks
pathways.show.all <- seurat.cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(seurat.cellchat, signaling = pathways.show.all[i])
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/LR_contribution_combined/L-R_contribution_combined_", pathways.show.all[i], ".pdf"), 
         plot = gg, 
         width = 13, 
         height = 7, 
         units = 'in', 
         dpi = 600)
}





## Generate violin plots for comparing Long Covid and Control for ligand-receptor pairs
# Extract all pathways
pathways.show.all <- seurat.cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  # Plot gene expression for each pathway in a violin plot, split by "disease"
  gg <- plotGeneExpression(seurat.cellchat, signaling = pathways.show.all[i], split.by = "disease", type = "violin")

  # Save each violin plot to a png file
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison/violin_LR_comparison_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}



# Plot gene expression for the TGFβ pathway in a violin plot, split by "disease"
gg <- plotGeneExpression(seurat.cellchat, signaling = "TGFb", split.by = "disease", angle.x = 45, type = "violin")

# Save the TGFβ pathway violin plot to a PNG file
ggsave(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison/violin_LR_comparison_TGFb_expression_02.png",
       plot = gg, width = 12, height = 7, dpi = 300)




## Generate violin plots for comparing Long Covid and Control for ligand-receptor pairs
pathways.show.all <- seurat.cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  # Plot gene expression for each pathway in a violin plot, with customization options
  gg <- plotGeneExpression(seurat.cellchat, 
                           signaling = pathways.show.all[i], 
                           enriched.only = TRUE, 
                           type = "violin", 
                           split.by = "disease",
                           colors.ggplot = TRUE) +
                           ggtitle(paste("Expression Levels for Signaling Pathway Components:", pathways.show.all[i]))

  
  # Save each violin plot to a PNG file
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison_02/violin_LR_comparison_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}


# Define the pathways of interest
pathways.interest <- c("JAM", "BTLA", "FLT3")

for (i in 1:length(pathways.interest)) {
  # Plot gene expression for each selected pathway in a violin plot
  gg <- plotGeneExpression(seurat.cellchat, 
                           signaling = pathways.interest[i], 
                           enriched.only = TRUE, 
                           type = "violin", 
                           split.by = "disease",
                           colors.ggplot = TRUE) +
                           ggtitle(paste("Expression Levels for Signaling Pathway Components:", pathways.interest[i]))

  # Save each violin plot to a PNG file
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison_03/violin_LR_comparison_", pathways.interest[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}





# Compute and visualize the network centrality scores
seurat.cellchat <- netAnalysis_computeCentrality(seurat.cellchat, slot.name = "netP")

# Signaling role analysis (scatter plot)
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_scatter_combined.png", 
    width = 5,
    height = 5, 
    units = "in",
    res = 500)
netAnalysis_signalingRole_scatter(seurat.cellchat)
dev.off()

# Signaling role analysis (heatmap)
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_heatmap_combined.png", 
    width = 12,
    height = 12, 
    units = "in",
    res = 500)
# Change the signaling role analysis heatmap font size using the 'font.size' argument
netAnalysis_signalingRole_heatmap(seurat.cellchat, pattern = "outgoing", font.size = 5) +
netAnalysis_signalingRole_heatmap(seurat.cellchat, pattern = "incoming", font.size = 5)

dev.off()


# Identify outgoing communication patterns of secreting cells
selectK(seurat.cellchat, pattern = "outgoing", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_outgoing_combined.png")
nPatterns = 6
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/outgoing_patterns_combined.png", 
    width = 10,
    height = 15, 
    units = "in",
    res = 300)
seurat.cellchat <- identifyCommunicationPatterns(seurat.cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

# Identify incoming communication patterns of target cells
selectK(seurat.cellchat, pattern = "incoming", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_incoming_combined.png")
nPatterns = 4
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/incoming_patterns_combined.png", 
    width = 10,
    height = 15, 
    units = "in",
    res = 300)
seurat.cellchat <- identifyCommunicationPatterns(seurat.cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()


library(ggalluvial)

# Generate river plot of signaling patterns
outgoing.signaling <- as.data.frame(seurat.cellchat@netP[["pattern"]][["outgoing"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
outgoing.signaling$Pattern <- factor(outgoing.signaling$Pattern, levels = paste0("Pattern ", 1:6))

incoming.signaling <- as.data.frame(seurat.cellchat@netP[["pattern"]][["incoming"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
incoming.signaling$Pattern <- factor(incoming.signaling$Pattern, levels = paste0("Pattern ", 1:4))

png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/netAnalysis_river_combined.png", 
    width = 13,
    height = 15, 
    units = "in",
    res = 500)
netAnalysis_river(seurat.cellchat, pattern = "outgoing") +
netAnalysis_river(seurat.cellchat, pattern = "incoming")
dev.off()
