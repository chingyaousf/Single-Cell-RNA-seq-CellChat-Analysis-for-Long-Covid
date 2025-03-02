# Load required libraries
rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)

## top20_genes_all_cell_types_with_symbols, keep all column (ID, Description, setSize,	enrichmentScore, NES, pvalue,	p.adjust,	qvalue,	rank, genes, leading_edge, group)
# Define the directory for reading RDS files
rds_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04"

# Define the directory for saving the output file
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes"

# Define the output Excel file path
output_excel_file <- file.path(output_directory, "top20_genes_all_cell_types_with_symbols_04.xlsx")


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
  rds_file <- file.path(rds_directory, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds"))
  
  if (file.exists(rds_file)) {
    # Load the RDS file
    pathway_enrichment_comp <- readRDS(rds_file)
    
    # Extract the GSEA results
    if ("gsea_results" %in% names(pathway_enrichment_comp)) {
      pathway_enrichment_df <- pathway_enrichment_comp$gsea_results
      
      # Check if the core_enrichment column exists
      if ("core_enrichment" %in% colnames(pathway_enrichment_df)) {
        # Extract the top 20 pathways by absolute NES
        top20_pathways <- pathway_enrichment_df %>%
          arrange(desc(abs(NES))) %>%   # Sort by absolute value of NES in descending order
          slice_head(n = 20) %>%       # Take the top 20 pathways
          mutate(genes = sapply(strsplit(as.character(core_enrichment), "/"), function(x) {
            # Convert Entrez IDs to gene symbols
            entrez_to_symbol <- bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
            paste(entrez_to_symbol$SYMBOL, collapse = ",") # Combine symbols into a single string
          })) %>% 
          dplyr::select(ID, Description, setSize,	enrichmentScore, NES, pvalue,	p.adjust,	qvalue,	rank, genes, leading_edge, group) # Select relevant columns
        
        # Create a sheet for this cell type
        addWorksheet(wb, cell_type)
        writeData(wb, cell_type, top20_pathways)
        
        message(paste("Added sheet for cell type:", cell_type))
      } else {
        message(paste("No 'core_enrichment' column found for cell type:", cell_type))
      }
    } else {
      message(paste("No 'gsea_results' element found in the RDS file for cell type:", cell_type))
    }
  } else {
    message(paste("RDS file not found for cell type:", cell_type))
  }
}

# Save the Excel file
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message("Excel file created successfully at: ", output_excel_file)






# Load required libraries
rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(openxlsx)
library(org.Hs.eg.db)

## top20_genes_all_cell_types_with_symbols, only keep column (Description, genes, group)
# Define the directory for reading RDS files
rds_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04"

# Define the directory for saving the output file
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes"

# Define the output Excel file path
output_excel_file <- file.path(output_directory, "top20_genes_all_cell_types_with_symbols_04_02.xlsx")


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
  rds_file <- file.path(rds_directory, paste0("pathway_enrichment_reactome_FindMarkers_", cell_type, ".rds"))
  
  if (file.exists(rds_file)) {
    # Load the RDS file
    pathway_enrichment_comp <- readRDS(rds_file)
    
    # Extract the GSEA results
    if ("gsea_results" %in% names(pathway_enrichment_comp)) {
      pathway_enrichment_df <- pathway_enrichment_comp$gsea_results
      
      # Check if the core_enrichment column exists
      if ("core_enrichment" %in% colnames(pathway_enrichment_df)) {
        # Extract the top 20 pathways by absolute NES
        top20_pathways <- pathway_enrichment_df %>%
          arrange(desc(abs(NES))) %>%   # Sort by absolute value of NES in descending order
          slice_head(n = 20) %>%       # Take the top 20 pathways
          mutate(genes = sapply(strsplit(as.character(core_enrichment), "/"), function(x) {
            # Convert Entrez IDs to gene symbols
            entrez_to_symbol <- bitr(x, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
            paste(entrez_to_symbol$SYMBOL, collapse = ",") # Combine symbols into a single string
          })) %>% 
          dplyr::select(Description, genes, group) # Select relevant columns
        
        # Create a sheet for this cell type
        addWorksheet(wb, cell_type)
        writeData(wb, cell_type, top20_pathways)
        
        message(paste("Added sheet for cell type:", cell_type))
      } else {
        message(paste("No 'core_enrichment' column found for cell type:", cell_type))
      }
    } else {
      message(paste("No 'gsea_results' element found in the RDS file for cell type:", cell_type))
    }
  } else {
    message(paste("RDS file not found for cell type:", cell_type))
  }
}

# Save the Excel file
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message("Excel file created successfully at: ", output_excel_file)










rds_file <- readRDS(file="/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/pathway_enrichment_reactome_FindMarkers_B_exhausted.rds")
str(rds_file)

rds_file <- readRDS(file="/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/GSEA_analysis/Reactome_GSEA_exclude_MT_RPL_RPS_genes/Reactome_GSEA_exclude_MT_RPL_RPS_genes_04/pathway_enrichment_reactome_FindMarkers_CD4.Tfh.rds")
str(rds_file)

# Extract the GSEA results
gsea_results <- rds_file$gsea_results

# Convert to a data frame to check the first few rows
gsea_df <- as.data.frame(gsea_results)
head(gsea_df)
