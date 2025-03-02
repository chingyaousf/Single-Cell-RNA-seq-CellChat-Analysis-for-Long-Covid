# Load required libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)


## signaling_pathways_matched_Reactome
# File paths
signaling_pathways_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signaling_pathways_genes.xlsx"
reactome_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/top15_genes_longcovid_control_all_cell_types_02.xlsx"
output_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/signaling_pathways_matched_reactome.xlsx"

# Load signaling pathways file
signaling_pathways <- read_excel(signaling_pathways_file) %>%
  mutate(Genes = str_split(Genes, ",\\s*")) # Split genes into lists

# Get all sheet names from the Reactome file                                
sheet_names <- excel_sheets(reactome_file)

# Function to match genes
match_genes <- function(signaling_genes, reactome_genes) {
  intersect(unlist(signaling_genes), unlist(reactome_genes))
}

# Initialize a list to store results
results <- list()

# Loop through each sheet (cell type)
for (sheet in sheet_names) {
  # Load data from the current sheet and check column names
  reactome_data <- read_excel(reactome_file, sheet = sheet, col_names = TRUE)
  
  # Check column names
  print(paste("Columns in sheet", sheet, ":", paste(colnames(reactome_data), collapse = ", ")))
  
  # Use the correct column names for "Group", "Description", and "genes"
  reactome_data <- reactome_data %>%
    rename(Group = 1, Description = 2, genes = 3) %>% # Columns are numeric indices
    mutate(genes = str_split(genes, ",\\s*")) # Split genes into lists
  
  # Perform matching for each signaling pathway
  for (i in 1:nrow(signaling_pathways)) {
    pathway_name <- signaling_pathways$`Signaling Pathways`[i]
    signaling_genes <- signaling_pathways$Genes[[i]]
    
    matched <- reactome_data %>%
      rowwise() %>%
      mutate(Matched_Genes = list(match_genes(signaling_genes, genes))) %>%
      filter(length(Matched_Genes) > 0) %>%
      mutate(Signaling_Pathway = pathway_name, Cell_Type = sheet)
    
    # Append results
    if (nrow(matched) > 0) {
      results <- append(results, list(matched))
    }
  }
}

# Combine all results
final_results <- bind_rows(results)

# Save the results to an Excel file
write.xlsx(final_results, output_file)

cat("Results saved to:", output_file, "\n")



## signaling_pathways_matched_c2_cp
# Load required libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# File paths
signaling_pathways_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signaling_pathways_genes.xlsx"
c2_cp_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/c2_cp_analysis/top15_C2_CP_longcovid_control_all_cell_types_02.xlsx"
output_file_c2_cp <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/signaling_pathways_matched_C2_CP.xlsx"

# Load signaling pathways file
signaling_pathways <- read_excel(signaling_pathways_file) %>%
  mutate(Genes = str_split(Genes, ",\\s*")) # Split genes into lists

# Get all sheet names from the C2_CP file                                
sheet_names <- excel_sheets(c2_cp_file)

# Function to match genes
match_genes <- function(signaling_genes, pathway_genes) {
  intersect(unlist(signaling_genes), unlist(pathway_genes))
}

# Initialize a list to store results
results <- list()

# Loop through each sheet (cell type)
for (sheet in sheet_names) {
  # Load data from the current sheet and check column names
  c2_cp_data <- read_excel(c2_cp_file, sheet = sheet, col_names = TRUE)
  
  # Check column names for debugging
  print(paste("Columns in sheet", sheet, ":", paste(colnames(c2_cp_data), collapse = ", ")))
  
  # Use the correct column names for "Group", "Description", and "genes"
  c2_cp_data <- c2_cp_data %>%
    rename(Group = 1, Description = 2, genes = 3) %>% # Columns are numeric indices
    mutate(genes = str_split(genes, ",\\s*")) # Split genes into lists
  
  # Perform matching for each signaling pathway
  for (i in 1:nrow(signaling_pathways)) {
    pathway_name <- signaling_pathways$`Signaling Pathways`[i]
    signaling_genes <- signaling_pathways$Genes[[i]]
    
    matched <- c2_cp_data %>%
      rowwise() %>%
      mutate(Matched_Genes = list(match_genes(signaling_genes, genes))) %>%
      filter(length(Matched_Genes) > 0) %>%
      mutate(Signaling_Pathway = pathway_name, Cell_Type = sheet)
    
    # Append results
    if (nrow(matched) > 0) {
      results <- append(results, list(matched))
    }
  }
}

# Combine all results
final_results <- bind_rows(results)

# Save the results to an Excel file
write.xlsx(final_results, output_file_c2_cp)

cat("Results saved to:", output_file_c2_cp, "\n")

