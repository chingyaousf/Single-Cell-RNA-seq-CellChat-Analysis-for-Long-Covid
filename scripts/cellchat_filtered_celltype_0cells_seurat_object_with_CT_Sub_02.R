
# Load required packages
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)


# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_CT_Sub.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_CT_Sub.rds")

levels(longcovid_cellchat@idents)
levels(control_cellchat@idents)


# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Calculate the centrality scores for the signaling roles
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Create a list of CellChat objects for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)

# Calculate the number of links for each object in the list
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

# Control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link))

# Set the minimum and maximum link weights to control the dot size
weight.MinMax <- range(num.link)

# Generate the signaling role scatter plots for both Long Covid and Control
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg)

# Save the combined plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_scatter_comparison_CT_Sub.png", 
    width = 10,
    height = 5, 
    units = "in",
    res = 500)

print(combined_plot)
dev.off()



# Generate signaling changes scatter plots "B_cells", "Monocytes" cell type between Long Covid and Control
# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_CT_Sub.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_CT_Sub.rds")

levels(longcovid_cellchat@idents)
levels(control_cellchat@idents)

# Update the CellChat objects if required
longcovid_cellchat <- updateCellChat(longcovid_cellchat)
control_cellchat <- updateCellChat(control_cellchat)

group.new <- union(levels(control_cellchat@idents), levels(longcovid_cellchat@idents))

longcovid_cellchat <- liftCellChat(longcovid_cellchat, group.new)
control_cellchat <- liftCellChat(control_cellchat, group.new)

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")


object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


# head(cellchat_combined@idents)
slotNames(cellchat_combined)
str(cellchat_combined)
unique(cellchat_combined@meta$ident)
levels(cellchat_combined@meta$ident)


# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD4.Tfh", "Lymph_prolif", "gdT", "pDC")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub/signalingChanges_", cell_type, "_comparison_CT_Sub.png"),
      width = 12,
      height = 7, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}


# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD4.Tfh", "Lymph_prolif", "gdT", "pDC")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub/signalingChanges_", cell_type, "_comparison_CT_Sub_01.png"),
      width = 20,
      height = 15, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}



# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD4.Tfh", "Lymph_prolif", "gdT", "pDC")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub/signalingChanges_", cell_type, "_comparison_CT_Sub_02.png"),
      width = 100,
      height = 100, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}




# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("B_exhausted", "B_immature", "B_malignant", "B_naive", "B_non-switched_memory", "B_switched_memory")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub_B_cells/signalingChanges_", cell_type, "_comparison_CT_Sub.png"),
      width = 12,
      height = 7, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}


# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("B_exhausted", "B_immature", "B_malignant", "B_naive", "B_non-switched_memory", "B_switched_memory")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub_B_cells/signalingChanges_", cell_type, "_comparison_CT_Sub_01.png"),
      width = 20,
      height = 15, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}



# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("B_exhausted", "B_immature", "B_malignant", "B_naive", "B_non-switched_memory", "B_switched_memory")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison_CT_Sub_B_cells/signalingChanges_", cell_type, "_comparison_CT_Sub_02.png"),
      width = 100,
      height = 100, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}










cell_type <- "pDC"
gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
  ggtitle(paste("Signaling Changes for Cell Type:", cell_type))

print(gg_cell_type)





longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat_CT_Sub.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat_CT_Sub.rds")

# Update the CellChat objects if required
longcovid_cellchat <- updateCellChat(longcovid_cellchat)
control_cellchat <- updateCellChat(control_cellchat)

group.new <- union(levels(control_cellchat@idents), levels(longcovid_cellchat@idents))

longcovid_cellchat <- liftCellChat(longcovid_cellchat, group.new)
control_cellchat <- liftCellChat(control_cellchat, group.new)

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")


object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)



## Generate violin plots for comparing Long Covid and Control for ligand-receptor pairs

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined.rds")

# Inspect the Seurat object
head(seurat_object@meta.data)
seurat_object
unique(seurat_object@meta.data$CT_Sub)

# Set up CellChat database
CellChatDB <- CellChatDB.human

# Set cell identities to the "CT_Sub_combined" column and specify the default assay
Idents(seurat_object) <- seurat_object$CT_Sub
DefaultAssay(seurat_object) <- "RNA"

# Create CellChat object
seurat.cellchat <- createCellChat(object = seurat_object, group.by = "CT_Sub")

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
write.csv(communication.df, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/communication_interactions_CT_Sub.csv", row.names = FALSE)

# Infer the cell-cell communication at a signaling pathway level
seurat.cellchat <- computeCommunProbPathway(seurat.cellchat)
saveRDS(seurat.cellchat, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/seurat_cellchat_CT_Sub.RDS")

# Load the saved CellChat object
seurat.cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/seurat_cellchat_CT_Sub.RDS")

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
  ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison_CT_Sub/violin_LR_comparison_CT_Sub_", pathways.show.all[i], "_expression.png"),
         plot = gg, width = 12, height = 7, dpi = 300)
}



# Specify the pathway to plot, for example "IGF"
pathway_to_plot <- "IGF"

# Generate the violin plot for the selected pathway
gg <- plotGeneExpression(seurat.cellchat, 
                         signaling = pathway_to_plot, 
                         enriched.only = TRUE, 
                         type = "violin", 
                         split.by = "disease",
                         colors.ggplot = TRUE) +
                         theme(legend.position = "right")

# Save the violin plot to a PNG file with increased width and height
ggsave(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/violin_LR_comparison_CT_Sub/violin_LR_comparison_CT_Sub_", pathway_to_plot, "_expression_02.png"),
       plot = gg, width = 12, height = 7, dpi = 300)
