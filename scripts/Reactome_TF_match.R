# Load required libraries
library(readxl)
library(dplyr)
library(openxlsx)


## Reactome_pathways_matched_TF
# Define input file paths
reactome_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/top15_genes_longcovid_control_all_cell_types_02.xlsx"
tf_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_filtered.xlsx"

# Define output file path
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes"
output_file_path <- file.path(output_directory, "Reactome_pathways_matched_TF.xlsx")

# Read the sheet names
reactome_sheets <- excel_sheets(reactome_file_path)
tf_sheets <- c("LongCovid Top 20 TFs", "Control Top 20 TFs")

# Function to calculate Jaccard Index
calculate_jaccard <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  if (union_len == 0) return(0)
  return(intersect_len / union_len)
}

# Create an Excel workbook for saving results
wb <- createWorkbook()

# Loop over Reactome sheets
for (cell_type in reactome_sheets) {
  # Load Reactome gene data
  reactome_df <- read_excel(reactome_file_path, sheet = cell_type)
  reactome_df$genes <- strsplit(reactome_df$genes, ",")  # Split gene lists into vectors
  reactome_df$genes <- lapply(reactome_df$genes, function(x) toupper(trimws(x)))  # Clean gene names
  
  # Loop over TF sheets
  for (tf_sheet in tf_sheets) {
    message(paste("Comparing Reactome sheet:", cell_type, "with TF sheet:", tf_sheet))
    
    # Load TF-target data
    tf_df <- read_excel(tf_file_path, sheet = tf_sheet)
    tf_df$target_genes <- strsplit(tf_df$target_genes, ",")
    tf_df$target_genes <- lapply(tf_df$target_genes, function(x) toupper(trimws(x)))  # Clean TF gene names
    
    # Initialize a data frame to store matches
    matches <- data.frame(
      Reactome_Pathway = character(),
      TF_Source = character(),
      Group = character(),
      Jaccard_Index = numeric(),
      Overlapping_Genes = character(),
      Reactome_Genes = character(),
      TF_Target_Genes = character(),
      stringsAsFactors = FALSE
    )
    
    # Compare Reactome pathways with TF target genes
    for (i in 1:nrow(reactome_df)) {
      for (j in 1:nrow(tf_df)) {
        # Calculate overlap and Jaccard Index
        overlap_genes <- intersect(reactome_df$genes[[i]], tf_df$target_genes[[j]])
        jaccard_index <- length(overlap_genes) / length(union(reactome_df$genes[[i]], tf_df$target_genes[[j]]))
        
        if (jaccard_index > 0.05) {  # Threshold for similarity
          matches <- bind_rows(matches, data.frame(
            Reactome_Pathway = reactome_df$Description[i],  # Assuming "Description" column exists
            TF_Source = tf_df$source[j],  # Assuming "source" column exists
            Group = tf_sheet,
            Jaccard_Index = jaccard_index,
            Overlapping_Genes = paste(overlap_genes, collapse = ", "),  # Store overlapping genes
            Reactome_Genes = paste(reactome_df$genes[[i]], collapse = ", "),
            TF_Target_Genes = paste(tf_df$target_genes[[j]], collapse = ", ")
          ))
        }
      }
    }
    
    # Abbreviate sheet name and save results
    if (nrow(matches) > 0) {
      abbreviated_tf_sheet <- gsub("LongCovid", "LC", gsub("Control", "Ctrl", tf_sheet))
      sheet_name <- substr(paste(cell_type, gsub(" ", "_", abbreviated_tf_sheet), sep = "_"), 1, 31)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, matches)
      message(paste("Added results for:", sheet_name))
    } else {
      message(paste("No matches found for Reactome sheet:", cell_type, "and TF sheet:", tf_sheet))
    }
  }
}

# Save the workbook
saveWorkbook(wb, output_file_path, overwrite = TRUE)
message("Results saved to: ", output_file_path)




## Reactome_pathways_matched_TF_diff
# Updated input paths
reactome_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/top15_genes_longcovid_control_all_cell_types_02.xlsx"
tf_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_diff.xlsx"
output_directory <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes"
output_file_path <- file.path(output_directory, "Reactome_pathways_matched_TF_diff.xlsx")

# Function to calculate Jaccard Index
calculate_jaccard <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  if (union_len == 0) return(0)
  return(intersect_len / union_len)
}

# Create an Excel workbook for saving results
wb <- createWorkbook()

# Read Reactome sheet names
reactome_sheets <- excel_sheets(reactome_file_path)

# Read TF-target data (only one sheet "Top 20 TF Differences")
tf_df <- read_excel(tf_file_path, sheet = "Top 20 TF Differences")
tf_df$target_genes <- strsplit(tf_df$target_genes, ",")
tf_df$target_genes <- lapply(tf_df$target_genes, function(x) toupper(trimws(x)))  # Clean gene names

# Loop over Reactome sheets
for (cell_type in reactome_sheets) {
  # Load Reactome gene data
  reactome_df <- read_excel(reactome_file_path, sheet = cell_type)
  reactome_df$genes <- strsplit(reactome_df$genes, ",")  # Split gene lists into vectors
  reactome_df$genes <- lapply(reactome_df$genes, function(x) toupper(trimws(x)))  # Clean gene names
  
  # Initialize a data frame to store matches
  matches <- data.frame(
    Reactome_Pathway = character(),
    TF_Source = character(),
    Jaccard_Index = numeric(),
    Overlapping_Genes = character(),
    Reactome_Genes = character(),
    TF_Target_Genes = character(),
    stringsAsFactors = FALSE
  )
  
  # Compare Reactome pathways with TF target genes
  for (i in 1:nrow(reactome_df)) {
    for (j in 1:nrow(tf_df)) {
      # Calculate overlap and Jaccard Index
      overlap_genes <- intersect(reactome_df$genes[[i]], tf_df$target_genes[[j]])
      jaccard_index <- length(overlap_genes) / length(union(reactome_df$genes[[i]], tf_df$target_genes[[j]]))
      
      if (jaccard_index > 0.02) {  # Threshold for similarity
        matches <- bind_rows(matches, data.frame(
          Reactome_Pathway = reactome_df$Description[i],  # Assuming "Description" column exists
          TF_Source = tf_df$source[j],                   # "source" column in TF data
          Jaccard_Index = jaccard_index,
          Overlapping_Genes = paste(overlap_genes, collapse = ", "),  # Overlapping genes
          Reactome_Genes = paste(reactome_df$genes[[i]], collapse = ", "),
          TF_Target_Genes = paste(tf_df$target_genes[[j]], collapse = ", ")
        ))
      }
    }
  }
  
  # Add results to the workbook if matches exist
  if (nrow(matches) > 0) {
    sheet_name <- substr(paste(cell_type, "TF_Overlap", sep = "_"), 1, 31)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, matches)
    message(paste("Added results for:", sheet_name))
  } else {
    message(paste("No matches found for Reactome sheet:", cell_type))
  }
}

# Save the workbook
saveWorkbook(wb, output_file_path, overwrite = TRUE)
message("Results saved to: ", output_file_path)
