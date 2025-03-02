# Load required libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)


## signaling_pathways_matched_TF
# File paths
signaling_pathways_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signaling_pathways_genes.xlsx"
tf_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_filtered.xlsx"
output_file_tf <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/signaling_pathways_matched_TF.xlsx"

# Load signaling pathways file
signaling_pathways <- read_excel(signaling_pathways_file) %>%
  mutate(Genes = str_split(Genes, ",\\s*")) # Split genes into lists

# Get all sheet names from the TF table file
sheet_names <- excel_sheets(tf_file)

# Function to match genes
match_genes <- function(signaling_genes, tf_target_genes) {
  intersect(unlist(signaling_genes), unlist(tf_target_genes))
}

# Initialize a list to store results
results <- list()

# Loop through each sheet (TF group: LongCovid and Control)
for (sheet in sheet_names) {
  # Load data from the current sheet
  tf_data <- read_excel(tf_file, sheet = sheet, col_names = TRUE) %>%
    mutate(target_genes = str_split(target_genes, ",\\s*")) # Split target_genes into lists
  
  # Perform matching for each signaling pathway
  for (i in 1:nrow(signaling_pathways)) {
    pathway_name <- signaling_pathways$`Signaling Pathways`[i]
    signaling_genes <- signaling_pathways$Genes[[i]]
    
    matched <- tf_data %>%
      rowwise() %>%
      mutate(Matched_Genes = list(match_genes(signaling_genes, target_genes))) %>%
      filter(length(Matched_Genes) > 0) %>%
      mutate(Signaling_Pathway = pathway_name, TF_Group = sheet)
    
    # Append results
    if (nrow(matched) > 0) {
      results <- append(results, list(matched))
    }
  }
}

# Combine all results
final_results <- bind_rows(results)

# Save the results to an Excel file
write.xlsx(final_results, output_file_tf)

cat("Results saved to:", output_file_tf, "\n")




## signaling_pathways_matched_TF_diff
# Load required libraries
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# File paths
signaling_pathways_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/signaling_pathways_genes.xlsx"
tf_diff_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/organized_top20_TF_targets_diff.xlsx"
output_file_tf_diff <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/signaling_pathways_matched_TF_diff.xlsx"

# Load signaling pathways file
signaling_pathways <- read_excel(signaling_pathways_file) %>%
  mutate(Genes = str_split(Genes, ",\\s*")) # Split genes into lists

# Load organized TF targets file
tf_diff_data <- read_excel(tf_diff_file, sheet = "Top 20 TF Differences") %>%
  mutate(target_genes = str_split(target_genes, ",\\s*")) # Split target_genes into lists

# Function to match genes
match_genes <- function(signaling_genes, tf_target_genes) {
  intersect(unlist(signaling_genes), unlist(tf_target_genes))
}

# Initialize results
results <- list()

# Perform matching for each signaling pathway
for (i in 1:nrow(signaling_pathways)) {
  pathway_name <- signaling_pathways$`Signaling Pathways`[i]
  signaling_genes <- signaling_pathways$Genes[[i]]
  
  matched <- tf_diff_data %>%
    rowwise() %>%
    mutate(Matched_Genes = list(match_genes(signaling_genes, target_genes))) %>%
    filter(length(Matched_Genes) > 0) %>%
    mutate(Signaling_Pathway = pathway_name)
  
  # Append results
  if (nrow(matched) > 0) {
    results <- append(results, list(matched))
  }
}

# Combine all results
final_results <- bind_rows(results)

# Save the results to an Excel file
write.xlsx(final_results, output_file_tf_diff)

cat("Results saved to:", output_file_tf_diff, "\n")
