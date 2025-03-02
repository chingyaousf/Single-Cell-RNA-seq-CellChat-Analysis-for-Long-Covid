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
cell_type <- "NK"

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
head(cell.markers, 50) 
head(cell.markers$gene) 


# Optionally save to a CSV file for inspection
write.csv(cell.markers, "NK_differential_genes.csv", row.names = FALSE)


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




