# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

## CD4.Tfh
## using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_03/"

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
  stop("No ranked genes available after resolving duplicates.")
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


# Save GSEA results as an RDS file
saveRDS(gsea_results, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.rds")))

# Save all pathways to a CSV file for further inspection
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  write.csv(as.data.frame(gsea_results), 
            file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.csv")), 
            row.names = FALSE)
}

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




## CD8.TE
## using FindMarkers instead of FindAllMarkers and using sortVar_med = median(sortVar, na.rm = TRUE)  # Use median to resolve ties
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD8.TE"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_03/"

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
  stop("No ranked genes available after resolving duplicates.")
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

# Save GSEA results as an RDS file
saveRDS(gsea_results, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.rds")))

# Save all pathways to a CSV file for further inspection
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  write.csv(as.data.frame(gsea_results), 
            file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_02.csv")), 
            row.names = FALSE)
}

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




## CD4.Tfh
## using # Update gene symbols for significant markers # significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(httr)  # Required for UpdateSymbolList function


# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD4.Tfh"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_03/"

# Subset for the current cell type using metadata filtering
current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)

# Extract all expressed genes as background
expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)

# Update gene symbols to the latest nomenclature
#expressed_genes <- UpdateSymbolList(symbols = expressed_genes, verbose = TRUE)

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
cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")

# Update gene symbols to the latest nomenclature
#cell.markers$gene <- UpdateSymbolList(symbols = cell.markers$gene, verbose = TRUE)

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
  stop("No ranked genes available after resolving duplicates.")
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

# Save GSEA results as an RDS file
saveRDS(gsea_results, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_03.rds")))

# Save all pathways to a CSV file for further inspection
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  write.csv(as.data.frame(gsea_results), 
            file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_03.csv")), 
            row.names = FALSE)
}

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
ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers_03.png")), 
       plot = plot, width = 8, height = 6, dpi = 300)

message(paste("Pathway enrichment completed for cell type:", cell_type))




## CD8.TE
## using # Update gene symbols for significant markers # significant_markers$gene <- UpdateSymbolList(symbols = significant_markers$gene, verbose = TRUE)
# Load the necessary libraries
library(Seurat)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(httr)  # Required for UpdateSymbolList function


# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell type to analyze
cell_type <- "CD8.TE"

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_03/"

# Subset for the current cell type using metadata filtering
current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)

# Extract all expressed genes as background
expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)

# Update gene symbols to the latest nomenclature
#expressed_genes <- UpdateSymbolList(symbols = expressed_genes, verbose = TRUE)

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
cell.markers <- tibble::rownames_to_column(cell.markers, var = "gene")

# Update gene symbols to the latest nomenclature
#cell.markers$gene <- UpdateSymbolList(symbols = cell.markers$gene, verbose = TRUE)

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
  stop("No ranked genes available after resolving duplicates.")
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

# Save GSEA results as an RDS file
saveRDS(gsea_results, file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_03.rds")))

# Save all pathways to a CSV file for further inspection
if (!is.null(gsea_results) && nrow(as.data.frame(gsea_results)) > 0) {
  write.csv(as.data.frame(gsea_results), 
            file = file.path(output_dir, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, "_03.csv")), 
            row.names = FALSE)
}

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
ggsave(filename = file.path(output_dir, paste0("Barplot_", cell_type, "_Reactome_FindMarkers_03.png")), 
       plot = plot, width = 8, height = 6, dpi = 300)

message(paste("Pathway enrichment completed for cell type:", cell_type))
