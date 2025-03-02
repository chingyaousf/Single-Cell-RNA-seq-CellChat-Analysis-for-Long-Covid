# Load required libraries
library(Seurat)  # For single-cell RNA-seq analysis
library(dplyr)   # For data manipulation

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define the cell types to analyze separately
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Output directory
output_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/DEG_27celltypes_exclude_MT_RPL_RPS_genes/"

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
  # Subset for the current cell type using metadata filtering
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 # max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Filter out RPL, RPS, and mitochondrial genes
  cell.markers <- cell.markers %>%
    filter(!grepl("^RPL", gene) & 
           !grepl("^RPS", gene) & 
           !grepl("^MT-", gene) & 
           !grepl("^MR-", gene))
  
  # Save the filtered cell.markers object for the current cell type
  saveRDS(cell.markers, file = paste0(output_dir, "FindAllMarkers_", cell_type, ".rds"))
  
  message(paste("Saved filtered markers for cell type:", cell_type))
}
