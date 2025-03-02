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
unique(seurat_object@meta.data$disease)


# Subset the Seurat object for LongCovid
longcovid_seurat <- subset(seurat_object, subset = disease == "LongCovid")

# Set up CellChat database
CellChatDB <- CellChatDB.human

# Set cell identities to the "CT_Sub_combined" column and specify the default assay
Idents(longcovid_seurat) <- longcovid_seurat$CT_Sub_combined
DefaultAssay(longcovid_seurat) <- "RNA"

# Create CellChat object for LongCovid
longcovid_cellchat <- createCellChat(object = longcovid_seurat, group.by = "CT_Sub_combined")

# Drop unused factor levels from cell identities
longcovid_cellchat@idents <- droplevels(longcovid_cellchat@idents)

# Check if the data matrix matches between Seurat and CellChat objects
summary(longcovid_cellchat@data@i == longcovid_seurat@assays[["RNA"]]@data@i)['FALSE'] # NA if no false

# Set the CellChatDB for signaling
longcovid_cellchat@DB <- CellChatDB

# Pre-processing
longcovid_cellchat <- subsetData(longcovid_cellchat)
longcovid_cellchat <- identifyOverExpressedGenes(longcovid_cellchat)
longcovid_cellchat <- identifyOverExpressedInteractions(longcovid_cellchat)

# Compute the communication probability
longcovid_cellchat <- computeCommunProb(longcovid_cellchat)

# Infer cellular communication network
longcovid_cellchat <- filterCommunication(longcovid_cellchat, min.cells = 10)
communication.df <- subsetCommunication(longcovid_cellchat)

# Save the Data Frame
write.csv(communication.df, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/communication_interactions_longcovid.csv", row.names = FALSE)

# Infer the cell-cell communication at a signaling pathway level
longcovid_cellchat <- computeCommunProbPathway(longcovid_cellchat)
saveRDS(longcovid_cellchat, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")

# Load the saved CellChat object
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")


# Calculate the aggregated cell-cell communication network
longcovid_cellchat <- aggregateNet(longcovid_cellchat)

# Adjust the interaction matrices and proceed with visualizations
matrix <- as.data.frame(longcovid_cellchat@net$count)
matrix[matrix == 0] <- 0.000000001  # Replace zeros with a very small value
matrix <- data.matrix(matrix)

# Save adjusted matrix if needed
# write.csv(matrix, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_matrix_longcovid.csv", row.names = FALSE)

# Visualize the number of interactions with weights/strengths
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_number_longcovid.png", 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(longcovid_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, edge.width.max = 10, weight.scale = TRUE, label.edge = TRUE, vertex.label.cex = 2, title.name = "interaction number Long Covid")
dev.off()

png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_strength_longcovid.png", 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(longcovid_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, edge.width.max = 10, weight.scale = TRUE, label.edge = TRUE, vertex.label.cex = 2, title.name = "interaction strength Long Covid")
dev.off()


# Save plots of all inferred networks
pathways.show.all <- longcovid_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(longcovid_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/LR_contribution_longcovid/LR_contribution_longcovid_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}




# Generate violin plots for each pathway to examine ligand-receptor gene expression
pathways.show.all <- longcovid_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  # Plot gene expression for each pathway in a violin plot and save to a variable
  gg <- plotGeneExpression(longcovid_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")

  # Save each violin plot to a png file
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_long_covid/violin_LR_long_covid_", pathways.show.all[i], "_expression_longcovid.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}




# Compute and visualize the network centrality scores
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")

# Signaling role analysis (scatter plot)
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_scatter_longcovid.png", 
    width = 5, height = 5, units = "in", res = 500)
netAnalysis_signalingRole_scatter(longcovid_cellchat, title = "Signaling Role - Long Covid")
dev.off()

# Signaling role analysis (heatmap)
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_heatmap_longcovid.png", 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(longcovid_cellchat, pattern = "outgoing", font.size = 5, title = "Long Covid") +
netAnalysis_signalingRole_heatmap(longcovid_cellchat, pattern = "incoming", font.size = 5 , title = "Long Covid")
dev.off()

# Identify outgoing and incoming communication patterns for LongCovid
selectK(longcovid_cellchat, pattern = "outgoing", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_outgoing_longcovid.png")
nPatterns = 6
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/outgoing_patterns_longcovid.png", 
    width = 10, height = 15, units = "in", res = 300)
longcovid_cellchat <- identifyCommunicationPatterns(longcovid_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

selectK(longcovid_cellchat, pattern = "incoming", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_incoming_longcovid.png")
nPatterns = 4
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/incoming_patterns_longcovid.png", 
    width = 10, height = 15, units = "in", res = 300)
longcovid_cellchat <- identifyCommunicationPatterns(longcovid_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

library(ggalluvial)

# Generate river plot of signaling patterns
outgoing.signaling <- as.data.frame(longcovid_cellchat@netP[["pattern"]][["outgoing"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
outgoing.signaling$Pattern <- factor(outgoing.signaling$Pattern, levels = paste0("Pattern ", 1:6))

incoming.signaling <- as.data.frame(longcovid_cellchat@netP[["pattern"]][["incoming"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
incoming.signaling$Pattern <- factor(incoming.signaling$Pattern, levels = paste0("Pattern ", 1:4))

png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/netAnalysis_river_longcovid.png", 
    width = 13, height = 15, units = "in", res = 500)
netAnalysis_river(longcovid_cellchat, pattern = "outgoing") +
netAnalysis_river(longcovid_cellchat, pattern = "incoming")
dev.off()






# Subset the Seurat object for Control
control_seurat <- subset(seurat_object, subset = disease == "Control")

# Set cell identities to the "CT_Sub_combined" column and specify the default assay
Idents(control_seurat) <- control_seurat$CT_Sub_combined
DefaultAssay(control_seurat) <- "RNA"

# Create CellChat object for Control
control_cellchat <- createCellChat(object = control_seurat, group.by = "CT_Sub_combined")

# Drop unused factor levels from cell identities
control_cellchat@idents <- droplevels(control_cellchat@idents)

# Check if the data matrix matches between Seurat and CellChat objects
summary(control_cellchat@data@i == control_seurat@assays[["RNA"]]@data@i)['FALSE'] # NA if no false

# Set the CellChatDB for signaling
control_cellchat@DB <- CellChatDB

# Pre-processing for Control
control_cellchat <- subsetData(control_cellchat)
control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)

# Compute the communication probability for Control
control_cellchat <- computeCommunProb(control_cellchat)

# Infer cellular communication network for Control
control_cellchat <- filterCommunication(control_cellchat, min.cells = 10)
communication.df_control <- subsetCommunication(control_cellchat)

# Save the Data Frame for Control
write.csv(communication.df_control, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/communication_interactions_control.csv", row.names = FALSE)


# Infer the cell-cell communication at a signaling pathway level for Control
control_cellchat <- computeCommunProbPathway(control_cellchat)
saveRDS(control_cellchat, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")


# Load the saved CellChat object
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")


# Calculate the aggregated cell-cell communication network for Control
control_cellchat <- aggregateNet(control_cellchat)

# Adjust the interaction matrices and proceed with visualizations
matrix_control <- as.data.frame(control_cellchat@net$count)
matrix_control[matrix_control == 0] <- 0.000000001  # Replace zeros with a very small value
matrix_control <- data.matrix(matrix_control)

# Save adjusted matrix if needed
# write.csv(matrix_control, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_matrix_control.csv", row.names = FALSE)

# Visualize the number of interactions with weights/strengths for Control
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_number_control.png", 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(control_cellchat@net$count, top = 0.1, arrow.width = 2, arrow.size = 0.5, edge.width.max = 10, weight.scale = TRUE, label.edge = TRUE, vertex.label.cex = 2, title.name = "interaction number Control")
dev.off()

png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_strength_control.png", 
    width = 10, height = 12, units = "in", res = 500)
netVisual_circle(control_cellchat@net$weight, top = 0.1, arrow.width = 2, arrow.size = 0.5, edge.width.max = 10, weight.scale = TRUE, label.edge = TRUE, vertex.label.cex = 2, title.name = "interaction strength Control")
dev.off()

# Save plots of all inferred networks for Control
pathways.show.all <- control_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution(control_cellchat, signaling = pathways.show.all[i])
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/LR_contribution_control/LR_contribution_control_", pathways.show.all[i], ".pdf"), 
         plot = gg, width = 13, height = 7, units = 'in', dpi = 600)
}



# Generate violin plots for each pathway to examine ligand-receptor gene expression
pathways.show.all <- control_cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  # Plot gene expression for each pathway in a violin plot and save to a variable
  gg <- plotGeneExpression(control_cellchat, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")

  # Save each violin plot to a png file
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_control/violin_LR_control_", pathways.show.all[i], "_expression_control.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}




# Compute and visualize the network centrality scores for Control
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Signaling role analysis (scatter plot) for Control
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_scatter_control.png", 
    width = 5, height = 5, units = "in", res = 500)
netAnalysis_signalingRole_scatter(control_cellchat, title = "Signaling Role - Control")
dev.off()

# Signaling role analysis (heatmap) for Control
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_heatmap_control.png", 
    width = 12, height = 12, units = "in", res = 500)
netAnalysis_signalingRole_heatmap(control_cellchat, pattern = "outgoing", font.size = 5, title = "Control") +
netAnalysis_signalingRole_heatmap(control_cellchat, pattern = "incoming", font.size = 5, title = "Control")
dev.off()

# Identify outgoing and incoming communication patterns for Control
selectK(control_cellchat, pattern = "outgoing", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_outgoing_control.png")
nPatterns = 6
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/outgoing_patterns_control.png", 
    width = 10, height = 15, units = "in", res = 300)
control_cellchat <- identifyCommunicationPatterns(control_cellchat, pattern = "outgoing", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()

selectK(control_cellchat, pattern = "incoming", nrun = 5)
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/selectK_incoming_control.png")
nPatterns = 4
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/incoming_patterns_control.png", 
    width = 10, height = 15, units = "in", res = 300)
control_cellchat <- identifyCommunicationPatterns(control_cellchat, pattern = "incoming", k = nPatterns, width = 3, height = 20, font.size = 10)
dev.off()


# Generate river plot of signaling patterns for Control
outgoing.signaling <- as.data.frame(control_cellchat@netP[["pattern"]][["outgoing"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
outgoing.signaling$Pattern <- factor(outgoing.signaling$Pattern, levels = paste0("Pattern ", 1:6))

incoming.signaling <- as.data.frame(control_cellchat@netP[["pattern"]][["incoming"]][["pattern"]][["signaling"]]) %>% group_by(Signaling) %>% slice_max(n = 1, order_by = Contribution)
incoming.signaling$Pattern <- factor(incoming.signaling$Pattern, levels = paste0("Pattern ", 1:4))

png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/netAnalysis_river_control.png", 
    width = 13, height = 15, units = "in", res = 500)

netAnalysis_river(control_cellchat, pattern = "outgoing") +
netAnalysis_river(control_cellchat, pattern = "incoming")
dev.off()







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

# Calculate the aggregated cell-cell communication network
seurat.cellchat <- aggregateNet(seurat.cellchat)

# Subset CellChat for LongCovid and Control
longcovid_cellchat <- subsetData(seurat.cellchat, subset.name = "disease", criteria = "LongCovid")
control_cellchat <- subsetData(seurat.cellchat, subset.name = "disease", criteria = "Control")

# Store both objects in a list for easier processing
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)

# Now, adjust the interaction matrices and proceed with visualizations.
for (i in 1:length(object.list)) {
  # Adjust the interaction count matrix (replacing zeros with a very small value)
  matrix <- as.data.frame(object.list[[i]]@net$count)
  matrix[matrix == 0] <- 0.000000001
  matrix <- data.matrix(matrix)
  
  # Save the adjusted matrix if needed
  #write.csv(matrix, file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_matrix_", names(object.list)[i], ".csv"), row.names = FALSE)
  
  # Interaction Numbers Visualization
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_number_", names(object.list)[i], ".png"), 
      width = 10, height = 10, units = "in", res = 500)
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.width.max = 10, 
                   arrow.width = 2, 
                   arrow.size = 0.5,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
  dev.off()
  
  # Interaction Strengths Visualization
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/interaction_strength_", names(object.list)[i], ".png"), 
      width = 10, height = 10, units = "in", res = 500)
  netVisual_circle(object.list[[i]]@net$weight, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.width.max = 10, 
                   arrow.width = 2, 
                   arrow.size = 0.5,
                   title.name = paste0("Interaction strength - ", names(object.list)[i]))
  dev.off()
}






























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
