# Load required libraries
library(readxl)
library(dplyr)
library(openxlsx)

# Define input file paths for Reactome and C2:CP data
reactome_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/top15_genes_longcovid_control_all_cell_types_02.xlsx"
c2_cp_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/c2_cp_analysis/top15_C2_CP_longcovid_control_all_cell_types_02.xlsx"

# Define the output directory
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes"
output_file_path <- file.path(output_directory, "matched_pathways_B_naive_02.xlsx")

# Load Reactome and C2:CP data for the "B_naive" sheet
reactome_df <- read_excel(reactome_file_path, sheet = "B_naive")
c2_cp_df <- read_excel(c2_cp_file_path, sheet = "B_naive")

# Split gene lists into R lists
reactome_df$genes <- strsplit(reactome_df$genes, ",")
c2_cp_df$genes <- strsplit(c2_cp_df$genes, ",")

# Function to calculate Jaccard Index
calculate_jaccard <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  if (union_len == 0) {
    return(0)
  }
  return(intersect_len / union_len)
}

# Compare Reactome and C2:CP pathways
matches <- data.frame(
  Reactome_Pathway = character(),
  C2_CP_Pathway = character(),
  Group = character(),
  Jaccard_Index = numeric(),
  Reactome_Genes = character(),
  C2_CP_Genes = character(),
  C2_CP_Adjusted_P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(reactome_df)) {
  for (j in 1:nrow(c2_cp_df)) {
    jaccard_index <- calculate_jaccard(reactome_df$genes[[i]], c2_cp_df$genes[[j]])
    if (jaccard_index > 0.2) {  # Threshold for similarity
      matches <- rbind(matches, data.frame(
        Reactome_Pathway = reactome_df$Description[i],
        C2_CP_Pathway = c2_cp_df$Description[j],
        Group = reactome_df$Group[i],  # LongCovid or Control
        Jaccard_Index = jaccard_index,
        Reactome_Genes = paste(reactome_df$genes[[i]], collapse = ", "),
        C2_CP_Genes = paste(c2_cp_df$genes[[j]], collapse = ", "),
        C2_CP_Adjusted_P_Value = ifelse("p.adjust" %in% names(c2_cp_df), c2_cp_df$p.adjust[j], NA)
      ))
    }
  }
}

# Save the results to an Excel file
write.xlsx(matches, output_file_path, rowNames = FALSE)
message("Matched pathways saved to: ", output_file_path)




# Load required libraries
library(readxl)
library(dplyr)
library(openxlsx)

# Define input file paths for Reactome and C2:CP data
reactome_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/top15_genes_longcovid_control_all_cell_types_02.xlsx"
c2_cp_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/c2_cp_analysis/top15_C2_CP_longcovid_control_all_cell_types_02.xlsx"

# Define the output directory
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes"
output_file_path <- file.path(output_directory, "matched_pathways_all_cell_types.xlsx")

# Read the sheet names from both Reactome and C2:CP files
reactome_sheets <- excel_sheets(reactome_file_path)
c2_cp_sheets <- excel_sheets(c2_cp_file_path)

# Find common cell types between Reactome and C2:CP
common_cell_types <- intersect(reactome_sheets, c2_cp_sheets)

# Create an Excel workbook to save results
wb <- createWorkbook()

# Function to calculate Jaccard Index
calculate_jaccard <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  if (union_len == 0) {
    return(0)
  }
  return(intersect_len / union_len)
}

# Loop over all common cell types
for (cell_type in common_cell_types) {
  message(paste("Processing cell type:", cell_type))
  
  # Load data for the current cell type
  reactome_df <- read_excel(reactome_file_path, sheet = cell_type)
  c2_cp_df <- read_excel(c2_cp_file_path, sheet = cell_type)
  
  # Split gene lists into R lists
  reactome_df$genes <- strsplit(reactome_df$genes, ",")
  c2_cp_df$genes <- strsplit(c2_cp_df$genes, ",")
  
  # Initialize a data frame to store matches
  matches <- data.frame(
    Reactome_Pathway = character(),
    C2_CP_Pathway = character(),
    Group = character(),
    Jaccard_Index = numeric(),
    Reactome_Genes = character(),
    C2_CP_Genes = character(),
    C2_CP_Adjusted_P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Compare Reactome and C2:CP pathways
  for (i in 1:nrow(reactome_df)) {
    for (j in 1:nrow(c2_cp_df)) {
      jaccard_index <- calculate_jaccard(reactome_df$genes[[i]], c2_cp_df$genes[[j]])
      if (jaccard_index > 0.2) {  # Threshold for similarity
        matches <- rbind(matches, data.frame(
          Reactome_Pathway = reactome_df$Description[i],
          C2_CP_Pathway = c2_cp_df$Description[j],
          Group = reactome_df$Group[i],  # LongCovid or Control
          Jaccard_Index = jaccard_index,
          Reactome_Genes = paste(reactome_df$genes[[i]], collapse = ", "),
          C2_CP_Genes = paste(c2_cp_df$genes[[j]], collapse = ", "),
          C2_CP_Adjusted_P_Value = ifelse("p.adjust" %in% names(c2_cp_df), c2_cp_df$p.adjust[j], NA)
        ))
      }
    }
  }
  
  # Add results to the Excel workbook as a sheet
  if (nrow(matches) > 0) {
    addWorksheet(wb, cell_type)
    writeData(wb, cell_type, matches)
    message(paste("Added results for cell type:", cell_type))
  } else {
    message(paste("No matches found for cell type:", cell_type))
  }
}

# Save the Excel workbook
saveWorkbook(wb, output_file_path, overwrite = TRUE)
message("Matched pathways saved to: ", output_file_path)
