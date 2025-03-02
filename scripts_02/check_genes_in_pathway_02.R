# Load required libraries
library(clusterProfiler)
library(dplyr)
library(openxlsx)

# Define the directory and output file paths
rds_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes"
output_excel_file <- file.path(rds_directory, "top15_genes_longcovid_control_all_cell_types.xlsx") # Updated file name

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
  # Modify the RDS file path to include the new directory
  rds_file <- file.path(rds_directory, paste0("pathway_enrichment_reactome_", cell_type, ".rds"))
  
  if (file.exists(rds_file)) {
    # Load the RDS file
    pathway_enrichment_comp <- readRDS(rds_file)
    
    # Convert the pathway enrichment result to a data frame
    pathway_enrichment_df <- as.data.frame(pathway_enrichment_comp)
    
    # Check if the Cluster column exists
    if ("Cluster" %in% colnames(pathway_enrichment_df)) {
      # Filter for the top 15 pathways for LongCovid and Control
      top15_longcovid <- pathway_enrichment_df %>%
        filter(Cluster == "LongCovid") %>%
        arrange(p.adjust) %>%
        slice_head(n = 15) %>%
        mutate(Group = "LongCovid")
      
      top15_control <- pathway_enrichment_df %>%
        filter(Cluster == "Control") %>%
        arrange(p.adjust) %>%
        slice_head(n = 15) %>%
        mutate(Group = "Control")
      
      # Combine LongCovid and Control pathways into one sheet
      combined_data <- bind_rows(top15_longcovid, top15_control) %>%
        mutate(genes = sapply(strsplit(geneID, "/"), function(x) paste(x, collapse = ","))) %>%
        dplyr::select(Group, Description, genes) # Use dplyr::select to avoid conflicts
      
      # Create a sheet for this cell type
      addWorksheet(wb, cell_type)
      writeData(wb, cell_type, combined_data)
      
      message(paste("Added sheet for cell type:", cell_type))
    } else {
      message(paste("No 'Cluster' column found for cell type:", cell_type))
    }
  } else {
    message(paste("RDS file not found for cell type:", cell_type))
  }
}

# Save the Excel file
saveWorkbook(wb, output_excel_file, overwrite = TRUE)
message("Excel file created successfully at: ", output_excel_file)
