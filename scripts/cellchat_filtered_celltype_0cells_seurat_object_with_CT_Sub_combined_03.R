
rm(list = ls())

# Load required packages
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggrepel)




# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

#Check Data Consistency: Ensure that both longcovid_cellchat and control_cellchat have sufficient overlap in signaling pathways and cell types.
signaling_longcovid <- longcovid_cellchat@netP$pathways
signaling_control <- control_cellchat@netP$pathways

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

# Generate the signaling role scatter plots for both Long Covid and Control
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg)

# Save the combined plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingRole_scatter_comparison.png", 
    width = 10,
    height = 5, 
    units = "in",
    res = 500)

print(combined_plot)
dev.off()




# Generate signaling changes scatter plots "B_cells", "Monocytes" cell type between Long Covid and Control
# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

levels(longcovid_cellchat@idents)
levels(control_cellchat@idents)

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Merge the two CellChat objects into one for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Set the cell types to compare
cell_types_to_compare <- c("B_cells", "Monocytes")

# Generate signaling changes scatter plots for the selected cell types
gg_list <- list()
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_list[[i]] <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) 
    #ggtitle(paste("Signaling Changes for Cell Type:", cell_type)) 
    
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg_list)

# Save the combined plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_Monocytes_comparison.png",
    width = 12,
    height = 6, 
    units = "in",
    res = 300)

print(combined_plot)
dev.off()


# Generate signaling changes scatter plots for all cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD8_T_cells", "CD4_T_cells", "B_cells", "NK_cells", "Monocytes", 
                           "Platelets", "gdT_cells", "DC_cells", "T_cells", "Lymphocytes", 
                           "NKT_cells", "Plasma_cells", "RBCs", "Treg_cells", "Stem_cells")

# Generate signaling changes scatter plots for the selected cell types
gg_list <- list()
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_list[[i]] <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type)) 
}

# Combine the plots for easy comparison
combined_plot <- patchwork::wrap_plots(plots = gg_list, ncol = 3)

# Save the combined plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_All_CellTypes_comparison.png",
    width = 30,
    height = 30, 
    units = "in",
    res = 300)

print(combined_plot)
dev.off()


# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD8_T_cells", "CD4_T_cells", "B_cells", "NK_cells", "Monocytes", 
                           "Platelets", "gdT_cells", "DC_cells", "T_cells", "Lymphocytes", 
                           "NKT_cells", "Plasma_cells", "RBCs", "Treg_cells", "Stem_cells")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_", cell_type, "_comparison.png"),
      width = 20,
      height = 15, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}



# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD8_T_cells", "CD4_T_cells", "B_cells", "NK_cells", "Monocytes", 
                           "Platelets", "gdT_cells", "DC_cells", "T_cells", "Lymphocytes", 
                           "NKT_cells", "Plasma_cells", "RBCs", "Treg_cells", "Stem_cells")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_", cell_type, "_comparison_01.png"),
      width = 12,
      height = 7, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}






# Generate signaling changes scatter plots for each cell type between Long Covid and Control
# Set the cell types to compare
cell_types_to_compare <- c("CD8_T_cells", "CD4_T_cells", "B_cells", "NK_cells", "Monocytes", 
                           "Platelets", "gdT_cells", "DC_cells", "T_cells", "Lymphocytes", 
                           "NKT_cells", "Plasma_cells", "RBCs", "Treg_cells", "Stem_cells")

# Generate and save signaling changes scatter plot for each selected cell type
for (i in 1:length(cell_types_to_compare)) {
  cell_type <- cell_types_to_compare[i]
  
  # Create scatter plot for signaling changes in the selected cell type
  gg_cell_type <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = cell_type) +
    ggtitle(paste("Signaling Changes for Cell Type:", cell_type))
  
  # Save the plot as a PNG file
  png(filename = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_", cell_type, "_comparison_02.png"),
      width = 100,
      height = 100, 
      units = "in",
      res = 300)
  
  print(gg_cell_type)
  dev.off()
}



# Generate signaling changes scatter plot for "B_cells" cell type between Long Covid and Control
# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

levels(longcovid_cellchat@idents)
levels(control_cellchat@idents)

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Merge the two CellChat objects into one for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

levels(cellchat_combined@idents)

head(cellchat_combined@idents)
slotNames(cellchat_combined)
str(cellchat_combined)
unique(cellchat_combined@meta$ident)



# Generate signaling changes scatter plot 
gg_b_cells <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = "B_cells")+
#geom_text_repel(max.overlaps = Inf)

# Save the plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_comparison_02.png",
    width = 100,
    height = 100, 
    units = "in",
    res = 300)

print(gg_b_cells)
dev.off()




# Create a new ggplot without the shape aesthetic
gg_b_cells <- ggplot(plot_data, aes(x = x, y = y, label = label, color = colour)) +
  geom_point() +
  geom_text_repel(max.overlaps = 50) +
  ggtitle("Signaling Changes for Cell Type: B_cells") +
  xlab("Differential outgoing interaction strength") +
  ylab("Differential incoming interaction strength") +
  theme_minimal()

# Save the updated plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_comparison_no_shape.png",
    width = 6,
    height = 6, 
    units = "in",
    res = 300)

print(gg_b_cells)
dev.off()



# Generate signaling changes scatter plot for "B_cells"
gg_b_cells <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = "B_cells")

# Extract the data from the ggplot object to modify it
plot_data <- ggplot_build(gg_b_cells)$data[[1]]

# Add labels to the plot data for `geom_text_repel()`
if ("name" %in% colnames(plot_data)) {
  plot_data$label <- plot_data$name
} else if ("label" %in% colnames(plot_data)) {
  plot_data$label <- plot_data$label
} else {
  plot_data$label <- rownames(plot_data)
}

# Create a new ggplot using the extracted data and adding `geom_text_repel()` with jitter
gg_b_cells_jittered <- ggplot(plot_data, aes(x = x, y = y, label = label, color = factor(colour))) +
  geom_point(position = position_jitter(width = 0.05, height = 0.05)) +
  geom_text_repel(max.overlaps = 100) +
  ggtitle("Signaling Changes for Cell Type: B_cells (Jittered)") +
  xlab("Differential outgoing interaction strength") +
  ylab("Differential incoming interaction strength") +
  theme_minimal() +
  scale_color_discrete(name = "Interaction Type")

# Save the updated plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_comparison_jittered.png",
    width = 6,
    height = 6, 
    units = "in",
    res = 300)
print(gg_b_cells_jittered)
dev.off()


# Generate the original signaling changes scatter plot
gg_b_cells_original <- netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = "B_cells")
  
# Extract data from the ggplot object to modify it for better labeling
plot_data <- ggplot_build(gg_b_cells_original)$data[[1]]

# Add labels for `geom_text_repel()`
plot_data$label <- if ("label" %in% colnames(plot_data)) {
  plot_data$label
} else if ("text" %in% colnames(plot_data)) {
  plot_data$text
} else {
  rownames(plot_data)
}

# Create a new ggplot with jitter and alpha for better visualization of overlapping points
gg_b_cells_updated <- ggplot(plot_data, aes(x = x, y = y, label = label, color = factor(colour), shape = factor(shape))) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.05, height = 0.05)) +  # Add jitter and alpha transparency
  geom_text_repel(max.overlaps = 100) +  # Add text labels with increased max.overlaps to avoid overlaps
  ggtitle("Signaling Changes for Cell Type: B_cells") +
  xlab("Differential outgoing interaction strength") +
  ylab("Differential incoming interaction strength")

# Save the updated plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_comparison_alpha.png",
    width = 6,
    height = 6, 
    units = "in",
    res = 300)

print(gg_b_cells_updated)
dev.off()



# Create a new ggplot with jitter, alpha, and manual settings for legends
gg_b_cells_updated <- ggplot(plot_data, aes(x = x, y = y, label = label, color = factor(colour), shape = factor(shape))) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.05, height = 0.05)) +  # Add jitter and alpha transparency
  geom_text_repel(max.overlaps = 100) +  # Add text labels with increased max.overlaps to avoid overlaps
  ggtitle("Signaling Changes for Cell Type: B_cells") +
  xlab("Differential outgoing interaction strength") +
  ylab("Differential incoming interaction strength") +
  theme_minimal() +
  scale_color_manual(name = "Condition", values = c("LongCovid" = "#F8766D", "Control" = "grey10")) +
  scale_shape_manual(name = "Interaction Type", values = c("Shared" = 21, "Incoming Specific" = 22, "Outgoing Specific" = 24))

# Save the updated plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signalingChanges_comparison/signalingChanges_B_cells_comparison_alpha_updated.png",
    width = 6,
    height = 6, 
    units = "in",
    res = 300)
print(gg_b_cells_updated)
dev.off()






## Compare the overall information flow of each signaling pathway or ligand-receptor pair
# Load the saved CellChat objects for Long Covid and Control
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

# Merge the two CellChat objects into one for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))

# Compare the overall information flow of each signaling pathway or ligand-receptor pair
# Rank the signaling pathways based on the information flow between LongCovid and Control
# This will generate stacked and unstacked bar plots
gg1 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE) +
  ggtitle("Relative Information Flow: Long Covid vs Control") +
  theme(plot.title = element_text(size = 18, face = "bold"))

gg2 <- rankNet(cellchat_combined, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE) +
  ggtitle("Information Flow: Long Covid vs Control") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Display the plots together for comparison
combined_plot <- gg1 + gg2

# Save the combined plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/informationFlow_comparison.png",
    width = 16,
    height = 8, 
    units = "in",
    res = 300)

print(combined_plot)
dev.off()




## Identify signaling groups based on their functional similarity between Long Covid and Control
reticulate::import("umap")
future::plan("sequential")

# Load the saved CellChat objects for Long Covid and Control
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Merge the two CellChat objects into one for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))


# Identify signaling groups based on their functional similarity
ptm <- Sys.time()

# Compute network similarity based on functional similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "functional")
#> Compute signaling network similarity for datasets LongCovid and Control

# Perform manifold learning of the signaling networks
cellchat_combined <- netEmbedding(cellchat_combined, type = "functional")

# Manifold learning of the signaling networks for datasets LongCovid and Control

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "functional")


# Classification learning of the signaling networks for datasets LongCovid and Control

# Visualization in 2D-space
gg_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "functional", label.size = 3.5) +
  ggtitle("Functional Similarity between Signaling Networks: Long Covid vs Control") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the visualization as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/functional_similarity_LongCovid_Control.png",
    width = 10,
    height = 8, 
    units = "in",
    res = 300)

print(gg_similarity)
dev.off()

# Calculate execution time
execution_time <- Sys.time() - ptm
print(paste("Execution time:", execution_time))


#Check Data Consistency: Ensure that both longcovid_cellchat and control_cellchat have sufficient overlap in signaling pathways and cell types.
signaling_longcovid <- longcovid_cellchat@netP$pathways
signaling_control <- control_cellchat@netP$pathways

intersected_pathways <- intersect(signaling_longcovid, signaling_control)
print(intersected_pathways)


# Compute and visualize the pathway "functional" distance in the learned joint manifold
distance_similarity <- rankSimilarity(cellchat_combined, type = "functional")
#> Compute the distance of signaling networks between datasets LongCovid and Control

# Add a title to the distance similarity plot
distance_similarity <- distance_similarity +
  ggtitle("Distance Functional Similarity between Signaling Networks: Long Covid vs Control") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  # Increase x-axis label size
    axis.title.y = element_text(size = 30)   # Increase y-axis label size
  )

# Save the distance similarity plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/distance_functional_similarity_LongCovid_Control.png",
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()


## Identify signaling groups based on structure similarity

# Load the saved CellChat objects for Long Covid and Control
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

# Aggregate communication network for each dataset
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Compute centrality for each dataset
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Merge the two CellChat objects into one for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list))


# Compute network similarity based on structural similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "structural")
#> Compute signaling network similarity for datasets LongCovid and Control

# Perform manifold learning of the signaling networks based on structural similarity
cellchat_combined <- netEmbedding(cellchat_combined, type = "structural")
#> Manifold learning of the signaling networks for datasets LongCovid and Control

# Perform classification learning of the signaling networks
cellchat_combined <- netClustering(cellchat_combined, type = "structural")
#> Classification learning of the signaling networks for datasets LongCovid and Control

# Visualization in 2D-space
gg_structural_similarity <- netVisual_embeddingPairwise(cellchat_combined, type = "structural", label.size = 3.5) +
  ggtitle("Structural Similarity between Signaling Networks: Long Covid vs Control") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/structural_similarity_LongCovid_Control.png",
    width = 10,
    height = 8, 
    units = "in",
    res = 300)
print(gg_structural_similarity)
dev.off()

# Visualization with zoom-in
gg_structural_similarity_zoom <- netVisual_embeddingPairwiseZoomIn(cellchat_combined, type = "structural", nCol = 2) +
  ggtitle("Zoomed-in Structural Similarity between Signaling Networks: Long Covid vs Control") +
  theme(plot.title = element_text(size = 18, face = "bold"))

# Save the zoomed-in plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/structural_similarity_LongCovid_Control_zoom.png",
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(gg_structural_similarity_zoom)
dev.off()

# Compute and visualize the pathway structural distance in the learned joint manifold
distance_similarity <- rankSimilarity(cellchat_combined, type = "structural")
#> Compute the distance of signaling networks between datasets LongCovid and Control

# Add a title to the distance similarity plot
distance_similarity <- distance_similarity +
  ggtitle("Distance Structural Similarity between Signaling Networks: Long Covid vs Control") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 25),  # Increase x-axis label size
    axis.title.y = element_text(size = 30)   # Increase y-axis label size
  )

# Save the distance similarity plot as a PNG file
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/distance_structural_similarity_LongCovid_Control.png",
    width = 12,
    height = 8, 
    units = "in",
    res = 300)
print(distance_similarity)
dev.off()







## # Generate heatmaps for comparing incoming signaling patterns between Long Covid and Control
# Load required libraries
library(ComplexHeatmap)

# Load the saved CellChat objects
longcovid_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/longcovid_cellchat.rds")
control_cellchat <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/control_cellchat.rds")

# Prepare the CellChat objects by aggregating the communication networks
longcovid_cellchat <- aggregateNet(longcovid_cellchat)
control_cellchat <- aggregateNet(control_cellchat)

# Calculate the centrality scores for the signaling roles
longcovid_cellchat <- netAnalysis_computeCentrality(longcovid_cellchat, slot.name = "netP")
control_cellchat <- netAnalysis_computeCentrality(control_cellchat, slot.name = "netP")

# Create a list of CellChat objects for comparison
object.list <- list(LongCovid = longcovid_cellchat, Control = control_cellchat)

# Combine all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

# Define the index for the two CellChat objects
i = 1

# Combine all identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)

# Generate heatmaps for comparing incoming signaling patterns between Long Covid and Control
ht1_incoming <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 10)
ht2_incoming <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 10)

# Draw and save the combined incoming heatmap
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/heatmap_incoming_comparison.png",
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_incoming + ht2_incoming, ht_gap = unit(0.5, "cm"))
dev.off()

# Generate heatmaps for comparing outgoing signaling patterns between Long Covid and Control
ht1_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 10)
ht2_outgoing <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 10)

# Draw and save the combined outgoing heatmap
png(filename = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/heatmap_outgoing_comparison.png",
    width = 20,
    height = 20,
    units = "in",
    res = 500)

draw(ht1_outgoing + ht2_outgoing, ht_gap = unit(0.5, "cm"))
dev.off()





