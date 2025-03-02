# Load required libraries
library(clusterProfiler)
library(dplyr)
library(openxlsx)

# Define the directory and output file paths
rds_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes"
output_excel_file <- file.path(rds_directory, "top30_genes_all_cell_types_combined.xlsx") # Updated file name

# Create an Excel workbook
wb <- createWorkbook()

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant",
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM",
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE",
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT",
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Loop through each cell type and process the RDS files
for (cell_type in cell_types_to_analyze) {
  # Construct the RDS file path
  rds_file <- file.path(rds_directory, paste0("pathway_enrichment_reactome_", cell_type, ".rds"))
  
  if (file.exists(rds_file)) {
    # Load the RDS file
    pathway_enrichment_comp <- readRDS(rds_file)
    
    # Convert the pathway enrichment result to a data frame
    pathway_enrichment_df <- as.data.frame(pathway_enrichment_comp@result)
    
    # Check if the core_enrichment column exists
    if ("core_enrichment" %in% colnames(pathway_enrichment_df)) {
      # Extract the top 30 pathways
      top30_pathways <- pathway_enrichment_df %>%
        arrange(p.adjust) %>%   # Sort by adjusted p-value
        slice_head(n = 30) %>%  # Take the top 30 pathways
        mutate(genes = sapply(strsplit(as.character(core_enrichment), "/"), function(x) paste(x, collapse = ","))) %>% # Extract gene names
        dplyr::select(Description, NES, p.adjust, genes) # Select relevant columns
      
      # Create a sheet for this cell type
      addWorksheet(wb, cell_type)
      writeData(wb, cell_type, top30_pathways)
      
      message(paste("Added sheet for cell type:", cell_type))
    } else {
      message(paste("No 'core_enrichment' column found for cell type:", cell_type))
    }
  } else {
    message(paste("RDS file not found for cell type:", cell_type))
  }
}

# Save the Excel file
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message("Excel file created successfully at: ", output_excel_file)




# Load required libraries
library(clusterProfiler)
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)

# Define the directory and output file paths
rds_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes"
output_excel_file <- file.path(rds_directory, "top30_genes_all_cell_types_with_symbols.xlsx") # Updated file name

# Create an Excel workbook
wb <- createWorkbook()

# Define the cell types to analyze
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant",
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM",
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE",
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT",
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Loop through each cell type and process the RDS files
for (cell_type in cell_types_to_analyze) {
  # Construct the RDS file path
  rds_file <- file.path(rds_directory, paste0("pathway_enrichment_reactome_", cell_type, ".rds"))
  
  if (file.exists(rds_file)) {
    # Load the RDS file
    pathway_enrichment_comp <- readRDS(rds_file)
    
    # Convert the pathway enrichment result to a data frame
    pathway_enrichment_df <- as.data.frame(pathway_enrichment_comp@result)
    
    # Check if the core_enrichment column exists
    if ("core_enrichment" %in% colnames(pathway_enrichment_df)) {
      # Extract the top 30 pathways
      top30_pathways <- pathway_enrichment_df %>%
        arrange(p.adjust) %>%   # Sort by adjusted p-value
        slice_head(n = 30) %>%  # Take the top 30 pathways
        mutate(genes = sapply(strsplit(as.character(core_enrichment), "/"), function(x) {
          # Convert Entrez IDs to gene symbols
          entrez_to_symbol <- bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
          paste(entrez_to_symbol$SYMBOL, collapse = ",") # Combine symbols into a single string
        })) %>% 
        dplyr::select(Description, NES, p.adjust, genes) # Select relevant columns
      
      # Create a sheet for this cell type
      addWorksheet(wb, cell_type)
      writeData(wb, cell_type, top30_pathways)
      
      message(paste("Added sheet for cell type:", cell_type))
    } else {
      message(paste("No 'core_enrichment' column found for cell type:", cell_type))
    }
  } else {
    message(paste("RDS file not found for cell type:", cell_type))
  }
}

# Save the Excel file
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message("Excel file created successfully at: ", output_excel_file)








rds_file <- readRDS(file="/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_B_exhausted.rds")
str(rds_file)
