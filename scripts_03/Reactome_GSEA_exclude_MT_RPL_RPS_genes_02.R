
## the code including checking Filter out unwanted genes
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/"

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

# Run FindAllMarkers for the comparison of LongCovid vs Control
cell.markers <- FindAllMarkers(current_cell_type_cells, 
                               assay = "RNA", 
                               logfc.threshold = 0.25, 
                               max.cells.per.ident = 1000, 
                               only.pos = FALSE)


# View all differential genes
#print(cell.markers)
#head(cell.markers, 50) 
#head(cell.markers$gene) 


# Optionally save to a CSV file for inspection
write.csv(cell.markers, "CD4.Tfh_differential_genes.csv", row.names = FALSE)


# Check the total number of genes and count the "MR-" genes
total_genes <- nrow(cell.markers)
mr_genes <- cell.markers %>% filter(grepl("^MRPL|^MRPS", gene))

mr_genes[["gene"]]

# Print summary information
message(paste("Total number of genes before filtering:", total_genes))
message(paste("Number of 'MRPL,MPS' genes:", nrow(mr_genes)))

# If there are any "MRPL,MRPS" genes, display them
if (nrow(mr_genes) > 0) {
  print(mr_genes)
}

# show the counts of DE genes for each group (LongCovid vs Control)
table(cell.markers$cluster)

# A list of genes that are duplicated (appear in both groups)
duplicates <- cell.markers %>% 
  group_by(gene) %>% 
  filter(n() > 1)
print(duplicates)



# Filter out unwanted genes
filtered_genes <- cell.markers %>%
  filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
         !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
         !grepl("^MT-", gene) &  # Exclude mitochondrial genes
         !grepl("^MRPL|^MRPS", gene))   # Exclude mitochondrial ribosomal genes


# Check the total number of genes and count the "MR" genes after filtering
total_genes_after <- nrow(filtered_genes)
mr_genes_after <- filtered_genes %>% filter(grepl("^MRPL|^MRPS", gene))

# Print summary information after filtering
message(paste("Total number of genes after filtering:", total_genes_after))
message(paste("Number of 'MR' genes after filtering:", nrow(mr_genes_after)))



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
  stop(paste("No ranked genes available for cell type:", cell_type))
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
saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_03.rds"))

# Check if pathways are enriched
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  # Convert results to a data frame
  gsea_df <- as.data.frame(gsea_results)
  
  # Save all pathways to a CSV file for further inspection
  write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_03.csv"), row.names = FALSE)
  
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




## solved duplicates in ranked genes issue
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_02/"

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

# Run FindAllMarkers for the comparison of LongCovid vs Control
cell.markers <- FindAllMarkers(current_cell_type_cells, 
                               assay = "RNA", 
                               logfc.threshold = 0.25, 
                               max.cells.per.ident = 1000, 
                               only.pos = FALSE)

# Save the differential expression results for inspection
write.csv(cell.markers, file = file.path(output_dir, "CD4.Tfh_differential_genes.csv"), row.names = FALSE)


# Filter out unwanted genes
filtered_genes <- cell.markers %>%
  filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
         !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
         !grepl("^MT-", gene) &  # Exclude mitochondrial genes
         !grepl("^MRPL|^MRPS", gene))   # Exclude mitochondrial ribosomal genes

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

# Resolve duplicates by averaging positive and negative values and prioritizing positive in ties
resolved_significant_markers <- significant_markers %>%
  group_by(ENTREZID) %>%
  summarize(
    avg_positive = mean(sortVar[sortVar > 0], na.rm = TRUE),
    avg_negative = mean(sortVar[sortVar < 0], na.rm = TRUE),
    final_sortVar = case_when(
      !is.na(avg_positive) & !is.na(avg_negative) ~ 
        if_else(abs(avg_positive) >= abs(avg_negative), avg_positive, avg_negative),
      !is.na(avg_positive) ~ avg_positive,
      !is.na(avg_negative) ~ avg_negative,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(final_sortVar)) %>%
  ungroup()

# Prepare the ranked gene list
ranked_genes <- resolved_significant_markers$final_sortVar
names(ranked_genes) <- resolved_significant_markers$ENTREZID
ranked_genes <- sort(ranked_genes, decreasing = TRUE)



# Filter for significant markers and include the cluster column
significant_markers <- filtered_genes %>% 
  filter(p_val_adj < 0.05) %>% 
  dplyr::select(gene, sortVar, cluster)

# Convert gene symbols to Entrez IDs
gene_rank <- bitr(significant_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Merge the ranking variable and cluster with Entrez IDs
significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))

# Resolve duplicates by averaging positive and negative values, tracking the strongest signal cluster
resolved_significant_markers <- significant_markers %>%
  group_by(ENTREZID) %>%
  summarize(
    avg_positive = mean(sortVar[sortVar > 0], na.rm = TRUE),
    avg_negative = mean(sortVar[sortVar < 0], na.rm = TRUE),
    final_sortVar = case_when(
      !is.na(avg_positive) & !is.na(avg_negative) ~ 
        if_else(abs(avg_positive) >= abs(avg_negative), avg_positive, avg_negative),
      !is.na(avg_positive) ~ avg_positive,
      !is.na(avg_negative) ~ avg_negative,
      TRUE ~ NA_real_
    ),
    cluster = cluster[which.max(abs(sortVar))]  # Track the cluster with strongest signal
  ) %>%
  filter(!is.na(final_sortVar)) %>%
  ungroup()

# Adjust the sign of sortVar based on the cluster
resolved_significant_markers <- resolved_significant_markers %>%
  mutate(final_sortVar = if_else(cluster == "Control", -final_sortVar, final_sortVar))



# Ensure the ranked gene list is valid
if (is.null(ranked_genes) || length(ranked_genes) == 0) {
  stop("No ranked genes available after resolving duplicates.")
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
saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_03.rds"))

# Check if pathways are enriched
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  # Convert results to a data frame
  gsea_df <- as.data.frame(gsea_results)
  
  # Save all pathways to a CSV file for further inspection
  write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_04.csv"), row.names = FALSE)
  
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






## Adjust the sign of sortVar based on the cluster
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_02/"

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

# Run FindAllMarkers for the comparison of LongCovid vs Control
cell.markers <- FindAllMarkers(current_cell_type_cells, 
                               assay = "RNA", 
                               logfc.threshold = 0.25, 
                               max.cells.per.ident = 1000, 
                               only.pos = FALSE)

# Save the differential expression results for inspection
write.csv(cell.markers, file = file.path(output_dir, "CD4.Tfh_differential_genes.csv"), row.names = FALSE)


# Filter out unwanted genes
filtered_genes <- cell.markers %>%
  filter(!grepl("^RPL", gene) &  # Exclude ribosomal large subunit genes
         !grepl("^RPS", gene) &  # Exclude ribosomal small subunit genes
         !grepl("^MT-", gene) &  # Exclude mitochondrial genes
         !grepl("^MRPL|^MRPS", gene))   # Exclude mitochondrial ribosomal genes

# Add a pseudocount to p-values to avoid -Inf during log transformation
pseudocount <- 10e-300

# Create a ranking variable
filtered_genes$sortVar <- (-log10(filtered_genes$p_val_adj + pseudocount)) * sign(filtered_genes$avg_log2FC)


# Filter for significant markers and include the cluster column
significant_markers <- filtered_genes %>% 
  filter(p_val_adj < 0.05) %>% 
  dplyr::select(gene, sortVar, cluster)

# Convert gene symbols to Entrez IDs
gene_rank <- bitr(significant_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Merge the ranking variable and cluster with Entrez IDs
significant_markers <- left_join(significant_markers, gene_rank, by = c("gene" = "SYMBOL"))

# Resolve duplicates by averaging positive and negative values, tracking the strongest signal cluster
resolved_significant_markers <- significant_markers %>%
  group_by(ENTREZID) %>%
  summarize(
    avg_positive = mean(sortVar[sortVar > 0], na.rm = TRUE),
    avg_negative = mean(sortVar[sortVar < 0], na.rm = TRUE),
    final_sortVar = case_when(
      !is.na(avg_positive) & !is.na(avg_negative) ~ 
        if_else(abs(avg_positive) >= abs(avg_negative), avg_positive, avg_negative),
      !is.na(avg_positive) ~ avg_positive,
      !is.na(avg_negative) ~ avg_negative,
      TRUE ~ NA_real_
    ),
    cluster = cluster[which.max(abs(sortVar))]  # Track the cluster with strongest signal
  ) %>%
  filter(!is.na(final_sortVar)) %>%
  ungroup()

# Adjust the sign of sortVar based on the cluster
resolved_significant_markers <- resolved_significant_markers %>%
  mutate(final_sortVar = if_else(cluster == "Control", -final_sortVar, final_sortVar))



# Ensure the ranked gene list is valid
if (is.null(ranked_genes) || length(ranked_genes) == 0) {
  stop("No ranked genes available after resolving duplicates.")
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
saveRDS(gsea_results, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_04.rds"))

# Check if pathways are enriched
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  # Convert results to a data frame
  gsea_df <- as.data.frame(gsea_results)
  
  # Save all pathways to a CSV file for further inspection
  write.csv(gsea_df, file = paste0(output_dir, "pathway_enrichment_reactome_", cell_type, "_04.csv"), row.names = FALSE)
  
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
  ggsave(filename = paste0(output_dir, "Barplot_", cell_type, "_LCvsControl_Reactome_04.png"), 
         plot = plot, width = 8, height = 6, dpi = 300)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
} else {
  message(paste("No pathways enriched for cell type:", cell_type))
}



# checking # Adjust the sign of sortVar based on the cluster
table(resolved_significant_markers$cluster)

before <- resolved_significant_markers
resolved_significant_markers <- resolved_significant_markers %>%
  mutate(final_sortVar = if_else(cluster == "Control", -final_sortVar, final_sortVar))
after <- resolved_significant_markers

# Compare
all.equal(before$final_sortVar, after$final_sortVar)







# Find duplicate gene names
duplicates <- significant_markers %>% 
  filter(duplicated(ENTREZID) | duplicated(ENTREZID, fromLast = TRUE))

# Print the duplicate gene information
head(duplicates)

# Save the duplicates for further inspection
#write.csv(duplicates, file = "CD4.Tfh_Duplicate_Genes.csv", row.names = FALSE)



# Total rows in the duplicates dataframe
total_duplicate_rows <- nrow(duplicates)

# Total unique genes in the duplicates dataframe
total_unique_genes <- duplicates %>%
  distinct(gene) %>%
  nrow()

# Print the results
message(paste("Total rows with duplicates:", total_duplicate_rows))
message(paste("Total unique duplicated genes:", total_unique_genes))



# Total rows in ranked_genes (including unique and duplicate entries)
total_rows_ranked_genes <- length(ranked_genes)

# Total unique genes in the duplicates of ranked_genes
total_unique_genes_ranked <- length(unique_duplicated_genes_ranked)

# Print the results
message(paste("Total rows in ranked_genes:", total_rows_ranked_genes))
message(paste("Total unique duplicated genes in ranked_genes:", total_unique_genes_ranked))



# Updated Code with Additional Debugging for Ranked Genes
# Check the total rows and duplicates in ranked_genes
total_rows_ranked_genes <- length(ranked_genes)
total_unique_genes_ranked <- length(unique(names(ranked_genes)))

# Identify duplicates in the ranked_genes vector
duplicates_ranked_genes <- names(ranked_genes)[duplicated(names(ranked_genes))]

# Count the total number of duplicates
total_ranked_gene_duplicates <- length(duplicates_ranked_genes)

# Print summary of ranked_genes
message(paste("Total rows in ranked_genes:", total_rows_ranked_genes))
message(paste("Total unique genes in ranked_genes:", total_unique_genes_ranked))
message(paste("Total duplicates in ranked_genes:", total_ranked_gene_duplicates))

# Show a preview of duplicate Entrez IDs
if (total_ranked_gene_duplicates > 0) {
  message("Preview of duplicate Entrez IDs in ranked_genes:")
  print(head(duplicates_ranked_genes))
} else {
  message("No duplicates found in ranked_genes.")
}
