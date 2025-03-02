# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/"

# Loop through all cell types
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
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 # max.cells.per.ident = 1000, 
                                 only.pos = FALSE)
  
  # Filter out unwanted genes
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
           !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
           !grepl("^MT-", gene) &  # Exclude mitochondrial genes
           !grepl("^MR-", gene))   # Exclude mitochondrial ribosomal genes
  
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
  
  # Prepare the ranked gene list
  ranked_genes <- significant_markers$sortVar
  names(ranked_genes) <- significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  

  # Check if ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next  # Skip this cell type
  }


  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(ranked_genes, 
                             organism = "human", 
                             pvalueCutoff = 1, 
                             minGSSize = 10, 
                             maxGSSize = 200, 
                             eps = 0, 
                             nPermSimple = 100000, 
                             pAdjustMethod = "fdr", 
                             seed = TRUE)
  
  # Save GSEA results as an RDS file
  saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, ".rds"))
  
  # Check if pathways are enriched
  if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
    # Convert results to a data frame
    gsea_df <- as.data.frame(gsea_results)
    
    # Save all pathways to a CSV file for further inspection
    write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, ".csv"), row.names = FALSE)
    
    # Extract top pathways for plotting
    top_pathways <- gsea_df %>% 
      arrange(desc(abs(NES))) %>%   # Sort by absolute NES to capture both positive and negative scores
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
    ggsave(filename = paste0(output_dir, "Barplot_", cell_type, "_LCvsControl_Reactome.png"), 
           plot = plot, width = 8, height = 6, dpi = 300)
    
    message(paste("Pathway enrichment completed for cell type:", cell_type))
  } else {
    message(paste("No pathways enriched for cell type:", cell_type))
  }
}



# run individual cell type
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("CD4.Tfh", "CD8.TE")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/"

# Loop through all cell types
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
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 # max.cells.per.ident = 1000, 
                                 only.pos = FALSE)
  
  # Filter out unwanted genes
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
           !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
           !grepl("^MT-", gene) &  # Exclude mitochondrial genes
           !grepl("^MR-", gene))   # Exclude mitochondrial ribosomal genes
  
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
  
  # Prepare the ranked gene list
  ranked_genes <- significant_markers$sortVar
  names(ranked_genes) <- significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  

  # Check if ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next  # Skip this cell type
  }


  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(ranked_genes, 
                             organism = "human", 
                             pvalueCutoff = 1, 
                             minGSSize = 10, 
                             maxGSSize = 200, 
                             eps = 0, 
                             nPermSimple = 100000, 
                             pAdjustMethod = "fdr", 
                             seed = TRUE)
  
  # Save GSEA results as an RDS file
  saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, ".rds"))
  
  # Check if pathways are enriched
  if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
    # Convert results to a data frame
    gsea_df <- as.data.frame(gsea_results)
    
    # Save all pathways to a CSV file for further inspection
    write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, ".csv"), row.names = FALSE)
    
    # Extract top pathways for plotting
    top_pathways <- gsea_df %>% 
      arrange(desc(abs(NES))) %>%   # Sort by absolute NES to capture both positive and negative scores
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
    ggsave(filename = paste0(output_dir, "Barplot_", cell_type, "_LCvsControl_Reactome.png"), 
           plot = plot, width = 8, height = 6, dpi = 300)
    
    message(paste("Pathway enrichment completed for cell type:", cell_type))
  } else {
    message(paste("No pathways enriched for cell type:", cell_type))
  }
}

if (is.null(ranked_genes) || length(ranked_genes) == 0) {
  message(paste("No ranked genes available for cell type:", cell_type))
  next
}




# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/"

# Loop through all cell types
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
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = FALSE)
  
  # Filter out unwanted genes
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
           !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
           !grepl("^MT-", gene) &  # Exclude mitochondrial genes
           !grepl("^MR-", gene))   # Exclude mitochondrial ribosomal genes
  
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
  
  # Prepare the ranked gene list
  ranked_genes <- significant_markers$sortVar
  names(ranked_genes) <- significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  

  # Check if ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next  # Skip this cell type
  }


  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(ranked_genes, 
                             organism = "human", 
                             pvalueCutoff = 1, 
                             minGSSize = 10, 
                             maxGSSize = 200, 
                             eps = 0, 
                             nPermSimple = 100000, 
                             pAdjustMethod = "fdr", 
                             seed = TRUE)
  
  # Save GSEA results as an RDS file
  saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type,"_02",".rds"))
  
  # Check if pathways are enriched
  if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
    # Convert results to a data frame
    gsea_df <- as.data.frame(gsea_results)
    
    # Save all pathways to a CSV file for further inspection
    write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type,"_02",".csv"), row.names = FALSE)
    
    # Extract top pathways for plotting
    top_pathways <- gsea_df %>% 
      arrange(desc(abs(NES))) %>%   # Sort by absolute NES to capture both positive and negative scores
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
    ggsave(filename = paste0(output_dir, "Barplot_", cell_type, "_LCvsControl_Reactome_02.png"), 
           plot = plot, width = 8, height = 6, dpi = 300)
    
    message(paste("Pathway enrichment completed for cell type:", cell_type))
  } else {
    message(paste("No pathways enriched for cell type:", cell_type))
  }
}




# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze
cell_types_to_analyze <- c("CD4.Tfh","CD8.TE")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/"

# Loop through all cell types
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
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = FALSE)
  
  # Filter out unwanted genes
  filtered_genes <- cell.markers %>%
    filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
           !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
           !grepl("^MT-", gene) &  # Exclude mitochondrial genes
           !grepl("^MR-", gene))   # Exclude mitochondrial ribosomal genes
  
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
  
  # Prepare the ranked gene list
  ranked_genes <- significant_markers$sortVar
  names(ranked_genes) <- significant_markers$ENTREZID
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  

  # Check if ranked_genes is valid
  if (is.null(ranked_genes) || length(ranked_genes) == 0) {
    message(paste("No ranked genes available for cell type:", cell_type))
    next  # Skip this cell type
  }


  # Perform GSEA with Reactome pathways
  gsea_results <- gsePathway(ranked_genes, 
                             organism = "human", 
                             pvalueCutoff = 1, 
                             minGSSize = 10, 
                             maxGSSize = 200, 
                             eps = 0, 
                             nPermSimple = 100000, 
                             pAdjustMethod = "fdr", 
                             seed = TRUE)
  
  # Save GSEA results as an RDS file
  saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type,"_03",".rds"))
  
  # Check if pathways are enriched
  if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
    # Convert results to a data frame
    gsea_df <- as.data.frame(gsea_results)
    
    # Save all pathways to a CSV file for further inspection
    write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type,"_03",".csv"), row.names = FALSE)
    
    # Extract top pathways for plotting
    top_pathways <- gsea_df %>% 
      arrange(desc(abs(NES))) %>%   # Sort by absolute NES to capture both positive and negative scores
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
    ggsave(filename = paste0(output_dir, "Barplot_", cell_type, "_LCvsControl_Reactome_03.png"), 
           plot = plot, width = 8, height = 6, dpi = 300)
    
    message(paste("Pathway enrichment completed for cell type:", cell_type))
  } else {
    message(paste("No pathways enriched for cell type:", cell_type))
  }
}
