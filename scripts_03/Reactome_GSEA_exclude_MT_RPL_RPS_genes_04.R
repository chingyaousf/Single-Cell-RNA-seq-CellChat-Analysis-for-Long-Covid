# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)


## loop to all cell types, using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## save GSE output with group information (longcovid, control)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"


# Loop through each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Check if background genes were successfully converted
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Use FindMarkers for direct comparison of LongCovid vs Control
  cell.markers <- FindMarkers(
    current_cell_type_cells, 
    ident.1 = "LongCovid",  # Group of interest
    ident.2 = "Control",    # Reference group
    assay = "RNA", 
    logfc.threshold = 0.25, 
    only.pos = FALSE
  )
  
  # Add the gene column to the results
  cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")
  
  # Save the differential expression results for inspection
  write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MRPL|^MRPS", gene))
  
  # Add a pseudocount to p-values to avoid -Inf during log transformation
  pseudocount <- 10e-300
  
  # Create a ranking variable
  filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)
  
  # Filter for significant markers
  significant_markers <- filtered_genes %>% 
    filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene, sortVar)
  
  # Convert gene symbols to Entrez IDs
  gene_rank <- bitr(significant_markers$gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  
  # Merge the ranking variable with Entrez IDs
  significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))
  
  # Resolve duplicates by averaging sortVar for the same ENTREZID
  resolved_significant_markers <- significant_markers %>%
    group_by(ENTREZID) %>%
    summarize(
      sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
    ) %>%
    ungroup()
  
  # Prepare the ranked gene list
  ranked_genes <- resolved_significant_markers$sortVar_med
  names(ranked_genes) <- resolved_significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Ensure ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next
  }
  
  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(
    ranked_genes, 
    organism = "human", 
    pvalueCutoff = 1, 
    minGSSize = 10, 
    maxGSSize = 200, 
    eps = 0, 
    nPermSimple = 100000, 
    pAdjustMethod = "fdr", 
    seed = TRUE
  )
  
  # Save GSEA results with group information as metadata
  metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
  )

  # Add group-specific pathways to the GSEA results
  gsea_df <- as.data.frame(gsea_results)

  # Add group information based on NES
  gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

  # Save the updated data frame with group information to a CSV file
  write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".csv")), 
          row.names = FALSE)

  # Add the updated GSEA data frame to the metadata for the RDS file
  metadata$gsea_results <- gsea_df

  # Save the updated metadata as an RDS file
  saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds")))

  
  
  # Extract top pathways for plotting
  top_pathways <- as.data.frame(gsea_results) %>% 
    arrange(desc(abs(NES))) %>% 
    head(n = 20)
  
  # Create the bar plot
  plot <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", 
         title = paste("Top Enriched Pathways (", cell_type, ")")) +
    theme_minimal() +
    scale_fill_gradient(low = "#c11c84", high = "grey") +
    theme(panel.border = element_rect(colour = "gray70", fill = NA), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  # Save the plot with "Barplot" in the filename
  ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers.png")), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}




## Because no Treg pathway for avove so run only Treg, gdT testing, and foind no ranked genes in Treg 
## loop to (gdT, Treg) cell types, using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## save GSE output with group information (longcovid, control)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("Treg", "gdT")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"


# Loop through each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Check if background genes were successfully converted
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Use FindMarkers for direct comparison of LongCovid vs Control
  cell.markers <- FindMarkers(
    current_cell_type_cells, 
    ident.1 = "LongCovid",  # Group of interest
    ident.2 = "Control",    # Reference group
    assay = "RNA", 
    logfc.threshold = 0.25, 
    only.pos = FALSE
  )
  
  # Add the gene column to the results
  cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")
  
  # Save the differential expression results for inspection
  write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MRPL|^MRPS", gene))
  
  # Add a pseudocount to p-values to avoid -Inf during log transformation
  pseudocount <- 10e-300
  
  # Create a ranking variable
  filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)
  
  # Filter for significant markers
  significant_markers <- filtered_genes %>% 
    filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene, sortVar)
  
  # Convert gene symbols to Entrez IDs
  gene_rank <- bitr(significant_markers$gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  
  # Merge the ranking variable with Entrez IDs
  significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))
  
  # Resolve duplicates by averaging sortVar for the same ENTREZID
  resolved_significant_markers <- significant_markers %>%
    group_by(ENTREZID) %>%
    summarize(
      sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
    ) %>%
    ungroup()
  
  # Prepare the ranked gene list
  ranked_genes <- resolved_significant_markers$sortVar_med
  names(ranked_genes) <- resolved_significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Ensure ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next
  }
  
  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(
    ranked_genes, 
    organism = "human", 
    pvalueCutoff = 1, 
    minGSSize = 10, 
    maxGSSize = 200, 
    eps = 0, 
    nPermSimple = 100000, 
    pAdjustMethod = "fdr", 
    seed = TRUE
  )
  
  # Save GSEA results with group information as metadata
  metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
  )

  # Add group-specific pathways to the GSEA results
  gsea_df <- as.data.frame(gsea_results)

  # Add group information based on NES
  gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

  # Save the updated data frame with group information to a CSV file
  write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".csv")), 
          row.names = FALSE)

  # Add the updated GSEA data frame to the metadata for the RDS file
  metadata$gsea_results <- gsea_df

  # Save the updated metadata as an RDS file
  saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds")))

  
  
  # Extract top pathways for plotting
  top_pathways <- as.data.frame(gsea_results) %>% 
    arrange(desc(abs(NES))) %>% 
    head(n = 20)
  
  # Create the bar plot
  plot <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", 
         title = paste("Top Enriched Pathways (", cell_type, ")")) +
    theme_minimal() +
    scale_fill_gradient(low = "#c11c84", high = "grey") +
    theme(panel.border = element_rect(colour = "gray70", fill = NA), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  # Save the plot with "Barplot" in the filename
  ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers.png")), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}







# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)


## loop to all cell types, using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## using  # Update gene symbols for significant markers,  significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)
## save GSE output with group information (longcovid, control)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"


# Loop through each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Check if background genes were successfully converted
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Use FindMarkers for direct comparison of LongCovid vs Control
  cell.markers <- FindMarkers(
    current_cell_type_cells, 
    ident.1 = "LongCovid",  # Group of interest
    ident.2 = "Control",    # Reference group
    assay = "RNA", 
    logfc.threshold = 0.25, 
    only.pos = FALSE
  )
  
  # Add the gene column to the results
  cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")
  
  # Save the differential expression results for inspection
  #write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MRPL|^MRPS", gene))
  
  # Add a pseudocount to p-values to avoid -Inf during log transformation
  pseudocount <- 10e-300
  
  # Create a ranking variable
  filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)
  
  # Filter for significant markers
  significant_markers <- filtered_genes %>% 
    filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene, sortVar)


  # Update gene symbols for significant markers
  significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)
  
 
  # Convert gene symbols to Entrez IDs
  gene_rank <- bitr(significant_markers$gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  
  # Merge the ranking variable with Entrez IDs
  significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))
  
  # Resolve duplicates by averaging sortVar for the same ENTREZID
  resolved_significant_markers <- significant_markers %>%
    group_by(ENTREZID) %>%
    summarize(
      sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
    ) %>%
    ungroup()
  
  # Prepare the ranked gene list
  ranked_genes <- resolved_significant_markers$sortVar_med
  names(ranked_genes) <- resolved_significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Ensure ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next
  }
  
  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(
    ranked_genes, 
    organism = "human", 
    pvalueCutoff = 1, 
    minGSSize = 10, 
    maxGSSize = 200, 
    eps = 0, 
    nPermSimple = 100000, 
    pAdjustMethod = "fdr", 
    seed = TRUE
  )
  
 # Save GSEA results with group information as metadata
  metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
  )

  # Add group-specific pathways to the GSEA results
  gsea_df <- as.data.frame(gsea_results)

  # Add group information based on NES
  gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

  # Save the updated data frame with group information to a CSV file
  write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.csv")), 
          row.names = FALSE)

  # Add the updated GSEA data frame to the metadata for the RDS file
  metadata$gsea_results <- gsea_df

  # Save the updated metadata as an RDS file
  saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.rds")))

  
  
  # Extract top pathways for plotting
  top_pathways <- as.data.frame(gsea_results) %>% 
    arrange(desc(abs(NES))) %>% 
    head(n = 20)
  
  # Create the bar plot
  plot <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", 
         title = paste("Top Enriched Pathways (", cell_type, ")")) +
    theme_minimal() +
    scale_fill_gradient(low = "#c11c84", high = "grey") +
    theme(panel.border = element_rect(colour = "gray70", fill = NA), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  # Save the plot with "Barplot" in the filename
  ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers_02.png")), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}




## loop to only "CD16_mono", "CD14_mono", using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## Truncate or wrap pathway names for better visibility
## save GSE output with group information (longcovid, control)
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("CD16_mono", "CD14_mono")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"

# Loop through each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Check if background genes were successfully converted
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Use FindMarkers for direct comparison of LongCovid vs Control
  cell.markers <- FindMarkers(
    current_cell_type_cells, 
    ident.1 = "LongCovid",  # Group of interest
    ident.2 = "Control",    # Reference group
    assay = "RNA", 
    logfc.threshold = 0.25, 
    only.pos = FALSE
  )
  
  # Add the gene column to the results
  cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")
  
  # Save the differential expression results for inspection
  write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MRPL|^MRPS", gene))
  
  # Add a pseudocount to p-values to avoid -Inf during log transformation
  pseudocount <- 10e-300
  
  # Create a ranking variable
  filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)
  
  # Filter for significant markers
  significant_markers <- filtered_genes %>% 
    filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene, sortVar)
  
  # Convert gene symbols to Entrez IDs
  gene_rank <- bitr(significant_markers$gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  
  # Merge the ranking variable with Entrez IDs
  significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))
  
  # Resolve duplicates by averaging sortVar for the same ENTREZID
  resolved_significant_markers <- significant_markers %>%
    group_by(ENTREZID) %>%
    summarize(
      sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
    ) %>%
    ungroup()
  
  # Prepare the ranked gene list
  ranked_genes <- resolved_significant_markers$sortVar_med
  names(ranked_genes) <- resolved_significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Ensure ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next
  }
  
  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(
    ranked_genes, 
    organism = "human", 
    pvalueCutoff = 1, 
    minGSSize = 10, 
    maxGSSize = 200, 
    eps = 0, 
    nPermSimple = 100000, 
    pAdjustMethod = "fdr", 
    seed = TRUE
  )
  
  # Save GSEA results with group information as metadata
  metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
  )

  # Add group-specific pathways to the GSEA results
  gsea_df <- as.data.frame(gsea_results)

  # Add group information based on NES
  gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

  # Save the updated data frame with group information to a CSV file
  write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".csv")), 
          row.names = FALSE)

  # Add the updated GSEA data frame to the metadata for the RDS file
  metadata$gsea_results <- gsea_df

  # Save the updated metadata as an RDS file
  saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds")))

  
  # Extract top pathways for plotting
  top_pathways <- as.data.frame(gsea_results) %>% 
    arrange(desc(abs(NES))) %>% 
    head(n = 20)
  
  # Truncate or wrap pathway names for better visibility
  top_pathways$Description <- sapply(top_pathways$Description, function(x) {
    if (nchar(x) > 50) {
      paste(strwrap(x, width = 50), collapse = "\n") # Wrap long names
    } else {
      x
    }
  })
  
  # Create the bar plot
  plot <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", 
         title = paste("Top Enriched Pathways (", cell_type, ")")) +
    theme_minimal() +
    scale_fill_gradient(low = "#c11c84", high = "grey") +
    theme(panel.border = element_rect(colour = "gray70", fill = NA), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  # Save the plot with "Barplot" in the filename
  ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers.png")), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}





## loop to only "CD16_mono", "CD14_mono", using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## Update gene symbols for significant markers
## Truncate or wrap pathway names for better visibility
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("CD16_mono", "CD14_mono")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"

# Loop through each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
  
  # Check if background genes were successfully converted
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Use FindMarkers for direct comparison of LongCovid vs Control
  cell.markers <- FindMarkers(
    current_cell_type_cells, 
    ident.1 = "LongCovid",  # Group of interest
    ident.2 = "Control",    # Reference group
    assay = "RNA", 
    logfc.threshold = 0.25, 
    only.pos = FALSE
  )
  
  # Add the gene column to the results
  cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")
  
  # Save the differential expression results for inspection
  #write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MRPL|^MRPS", gene))
  
  # Add a pseudocount to p-values to avoid -Inf during log transformation
  pseudocount <- 10e-300
  
  # Create a ranking variable
  filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)
  
  # Filter for significant markers
  significant_markers <- filtered_genes %>% 
    filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene, sortVar)
  

  # Update gene symbols for significant markers
  significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)


  # Convert gene symbols to Entrez IDs
  gene_rank <- bitr(significant_markers$gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
  
  # Merge the ranking variable with Entrez IDs
  significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))
  
  # Resolve duplicates by averaging sortVar for the same ENTREZID
  resolved_significant_markers <- significant_markers %>%
    group_by(ENTREZID) %>%
    summarize(
      sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
    ) %>%
    ungroup()
  
  # Prepare the ranked gene list
  ranked_genes <- resolved_significant_markers$sortVar_med
  names(ranked_genes) <- resolved_significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Ensure ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next
  }
  
  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(
    ranked_genes, 
    organism = "human", 
    pvalueCutoff = 1, 
    minGSSize = 10, 
    maxGSSize = 200, 
    eps = 0, 
    nPermSimple = 100000, 
    pAdjustMethod = "fdr", 
    seed = TRUE
  )
  
  # Save GSEA results with group information as metadata
  metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
  )

  # Add group-specific pathways to the GSEA results
  gsea_df <- as.data.frame(gsea_results)

  # Add group information based on NES
  gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

  # Save the updated data frame with group information to a CSV file
  write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".csv")), 
          row.names = FALSE)

  # Add the updated GSEA data frame to the metadata for the RDS file
  metadata$gsea_results <- gsea_df

  # Save the updated metadata as an RDS file
  saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds")))

  
  # Extract top pathways for plotting
  top_pathways <- as.data.frame(gsea_results) %>% 
    arrange(desc(abs(NES))) %>% 
    head(n = 20)
  
  # Truncate or wrap pathway names for better visibility
  top_pathways$Description <- sapply(top_pathways$Description, function(x) {
    if (nchar(x) > 50) {
      paste(strwrap(x, width = 50), collapse = "\n") # Wrap long names
    } else {
      x
    }
  })
  
  # Create the bar plot
  plot <- ggplot(top_pathways, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score", 
         title = paste("Top Enriched Pathways (", cell_type, ")")) +
    theme_minimal() +
    scale_fill_gradient(low = "#c11c84", high = "grey") +
    theme(panel.border = element_rect(colour = "gray70", fill = NA), 
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
  
  # Save the plot with "Barplot" in the filename
  ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers_02.png")), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}



## CD4.Tfh cell types, using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## save GSE output with group information (longcovid, control)
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"

# Subset for the current cell type using metadata filtering
current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)

# Extract all expressed genes as background
expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)

# Convert expressed genes to Entrez IDs for background
background_entrez <- bitr(expressed_genes, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)

# Check if background genes were successfully converted
if (is.null(background_entrez) || nrow(background_entrez) == 0) {
  stop(paste("No background genes available for cell type:", cell_type))
}

# Set the identity to 'disease' for comparison
Idents(current_cell_type_cells) <- current_cell_type_cells$disease

# Use FindMarkers for direct comparison of LongCovid vs Control
cell.markers <- FindMarkers(
  current_cell_type_cells, 
  ident.1 = "LongCovid",  # Group of interest
  ident.2 = "Control",    # Reference group
  assay = "RNA", 
  logfc.threshold = 0.25, 
  only.pos = FALSE
)

# Add the gene column to the results
cell.markers <- rownames_to_column(cell.markers, var = "gene")

# Save the differential expression results for inspection
#write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)

# Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
filtered_genes <- cell.markers %>%
  filter(!grepl("^RPL", gene) & 
         !grepl("^RPS", gene) & 
         !grepl("^MT-", gene) & 
         !grepl("^MRPL|^MRPS", gene))

# Add a pseudocount to p-values to avoid -Inf during log transformation
pseudocount <- 10e-300

# Create a ranking variable
filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)

# Filter for significant markers
significant_markers <- filtered_genes %>% 
  filter(p_val_adj < 0.05) %>% 
  dplyr::select(gene, sortVar)

# Convert gene symbols to Entrez IDs
gene_rank <- bitr(significant_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Merge the ranking variable with Entrez IDs
significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))

# Resolve duplicates by averaging sortVar for the same ENTREZID
resolved_significant_markers <- significant_markers %>%
  group_by(ENTREZID) %>%
  summarize(
    sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
  ) %>%
  ungroup()

# Prepare the ranked gene list
ranked_genes <- resolved_significant_markers$sortVar_med
names(ranked_genes) <- resolved_significant_markers$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Ensure ranked_genes is valid
if (is.null(ranked_genes) || length(ranked_genes) == 0) {
  stop(paste("No ranked genes available for cell type:", cell_type))
}

# Perform GSEA with Reactome pathways
gsea_results <- gsePathway(
  ranked_genes, 
  organism = "human", 
  pvalueCutoff = 1, 
  minGSSize = 10, 
  maxGSSize = 200, 
  eps = 0, 
  nPermSimple = 100000, 
  pAdjustMethod = "fdr", 
  seed = TRUE
)

# Save GSEA results with group information as metadata
metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
)

# Add group-specific pathways to the GSEA results
gsea_df <- as.data.frame(gsea_results)

# Add group information based on NES
gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

# Save the updated data frame with group information to a CSV file
write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_05.csv")), 
          row.names = FALSE)

# Add the updated GSEA data frame to the metadata for the RDS file
metadata$gsea_results <- gsea_df

# Save the updated metadata as an RDS file
saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_05.rds")))

message(paste("Pathway enrichment completed for cell type:", cell_type))






## CD4.Tfh cell types, using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
## Update gene symbols for significant markers
## save GSE output with group information (longcovid, control)
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/"

# Subset for the current cell type using metadata filtering
current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)

# Extract all expressed genes as background
expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)

# Convert expressed genes to Entrez IDs for background
background_entrez <- bitr(expressed_genes, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)

# Check if background genes were successfully converted
if (is.null(background_entrez) || nrow(background_entrez) == 0) {
  stop(paste("No background genes available for cell type:", cell_type))
}

# Set the identity to 'disease' for comparison
Idents(current_cell_type_cells) <- current_cell_type_cells$disease

# Use FindMarkers for direct comparison of LongCovid vs Control
cell.markers <- FindMarkers(
  current_cell_type_cells, 
  ident.1 = "LongCovid",  # Group of interest
  ident.2 = "Control",    # Reference group
  assay = "RNA", 
  logfc.threshold = 0.25, 
  only.pos = FALSE
)

# Add the gene column to the results
cell.markers <- rownames_to_column(cell.markers, var = "gene")

# Save the differential expression results for inspection
#write.csv(cell.markers, file = file.path(output_dir, paste0(cell_type, "_differential_genes_FindMarkers.csv")), row.names = FALSE)

# Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
filtered_genes <- cell.markers %>%
  filter(!grepl("^RPL", gene) & 
         !grepl("^RPS", gene) & 
         !grepl("^MT-", gene) & 
         !grepl("^MRPL|^MRPS", gene))

# Add a pseudocount to p-values to avoid -Inf during log transformation
pseudocount <- 10e-300

# Create a ranking variable
filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)

# Filter for significant markers
significant_markers <- filtered_genes %>% 
  filter(p_val_adj < 0.05) %>% 
  dplyr::select(gene, sortVar)


# Update gene symbols for significant markers
  significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)



# Convert gene symbols to Entrez IDs
gene_rank <- bitr(significant_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Merge the ranking variable with Entrez IDs
significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))

# Resolve duplicates by averaging sortVar for the same ENTREZID
resolved_significant_markers <- significant_markers %>%
  group_by(ENTREZID) %>%
  summarize(
    sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
  ) %>%
  ungroup()

# Prepare the ranked gene list
ranked_genes <- resolved_significant_markers$sortVar_med
names(ranked_genes) <- resolved_significant_markers$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Ensure ranked_genes is valid
if (is.null(ranked_genes) || length(ranked_genes) == 0) {
  stop(paste("No ranked genes available for cell type:", cell_type))
}

# Perform GSEA with Reactome pathways
gsea_results <- gsePathway(
  ranked_genes, 
  organism = "human", 
  pvalueCutoff = 1, 
  minGSSize = 10, 
  maxGSSize = 200, 
  eps = 0, 
  nPermSimple = 100000, 
  pAdjustMethod = "fdr", 
  seed = TRUE
)

# Save GSEA results with group information as metadata
metadata <- list(
  gsea_results = gsea_results,  # GSEA result object
  comparison = c("LongCovid", "Control"),  # Group information
  cell_type = cell_type  # Cell type analyzed
)

# Add group-specific pathways to the GSEA results
gsea_df <- as.data.frame(gsea_results)

# Add group information based on NES
gsea_df$group <- ifelse(gsea_df$NES > 0, "LongCovid", "Control")

# Save the updated data frame with group information to a CSV file
write.csv(gsea_df, 
          file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_06.csv")), 
          row.names = FALSE)

# Add the updated GSEA data frame to the metadata for the RDS file
metadata$gsea_results <- gsea_df

# Save the updated metadata as an RDS file
saveRDS(metadata, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_06.rds")))

message(paste("Pathway enrichment completed for cell type:", cell_type))
