# Load the required libraries
library(dplyr)

# checking pathways, GenRatio, BgRatio
# File paths for Reactome and C2:CP results
reactome_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/comparasion_Reactome_c2cp/pathway_enrichment_reactome_B_exhausted.rds"
c2cp_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/comparasion_Reactome_c2cp/pathway_enrichment_c2_cp_B_exhausted.rds"

# Load Reactome and C2:CP results
reactome_results <- readRDS(reactome_file)  # Reactome enrichment results
c2cp_results <- readRDS(c2cp_file)         # C2:CP enrichment results

# Extract results as data frames from compareClusterResult objects
reactome_df <- reactome_results@compareClusterResult  # Use the slot for actual results
c2cp_df <- c2cp_results@compareClusterResult

# Check the structure of the extracted data frames
str(reactome_df)
str(c2cp_df)

# Filter significant pathways BEFORE and AFTER multiple testing correction
reactome_raw <- reactome_df %>% filter(pvalue <= 0.05)       # Filter based on raw p-value
reactome_adjusted <- reactome_df %>% filter(p.adjust <= 0.05) # Filter based on adjusted p-value

c2cp_raw <- c2cp_df %>% filter(pvalue <= 0.05)       # Filter based on raw p-value
c2cp_adjusted <- c2cp_df %>% filter(p.adjust <= 0.05) # Filter based on adjusted p-value

# Count the number of significant pathways
cat("Reactome - Significant pathways:\n")
cat("Raw p-value <= 0.05: ", nrow(reactome_raw), "\n")
cat("Adjusted p-value <= 0.05: ", nrow(reactome_adjusted), "\n")

cat("\nC2:CP - Significant pathways:\n")
cat("Raw p-value <= 0.05: ", nrow(c2cp_raw), "\n")
cat("Adjusted p-value <= 0.05: ", nrow(c2cp_adjusted), "\n")

# View top pathways for comparison
cat("\nTop Reactome Pathways:\n")
head(reactome_adjusted)

cat("\nTop C2:CP Pathways:\n")
head(c2cp_adjusted)




# Load the required libraries
library(dplyr)

## checking pathways, GenRatio, BgRatio
# File paths for Reactome, C2:CP, and C2 results
reactome_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/comparasion_Reactome_c2cp/pathway_enrichment_reactome_B_exhausted.rds"
c2cp_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/comparasion_Reactome_c2cp/pathway_enrichment_c2_cp_B_exhausted.rds"
c2_file <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/comparasion_Reactome_c2cp/pathway_enrichment_c2_B_exhausted.rds"

# Load Reactome, C2:CP, and C2 results
reactome_results <- readRDS(reactome_file)  # Reactome enrichment results
c2cp_results <- readRDS(c2cp_file)         # C2:CP enrichment results
c2_results <- readRDS(c2_file)            # C2 enrichment results

# Extract results as data frames from compareClusterResult objects
reactome_df <- reactome_results@compareClusterResult  # Reactome results
c2cp_df <- c2cp_results@compareClusterResult          # C2:CP results
c2_df <- c2_results@compareClusterResult             # C2 results

# Check the structure of the extracted data frames
str(reactome_df)
str(c2cp_df)
str(c2_df)

# Filter significant pathways BEFORE and AFTER multiple testing correction
reactome_raw <- reactome_df %>% filter(pvalue <= 0.05)       # Reactome raw p-value
reactome_adjusted <- reactome_df %>% filter(p.adjust <= 0.05) # Reactome adjusted p-value

c2cp_raw <- c2cp_df %>% filter(pvalue <= 0.05)       # C2:CP raw p-value
c2cp_adjusted <- c2cp_df %>% filter(p.adjust <= 0.05) # C2:CP adjusted p-value

c2_raw <- c2_df %>% filter(pvalue <= 0.05)           # C2 raw p-value
c2_adjusted <- c2_df %>% filter(p.adjust <= 0.05)    # C2 adjusted p-value

# Count the number of significant pathways
cat("Reactome - Significant pathways:\n")
cat("Raw p-value <= 0.05: ", nrow(reactome_raw), "\n")
cat("Adjusted p-value <= 0.05: ", nrow(reactome_adjusted), "\n")

cat("\nC2:CP - Significant pathways:\n")
cat("Raw p-value <= 0.05: ", nrow(c2cp_raw), "\n")
cat("Adjusted p-value <= 0.05: ", nrow(c2cp_adjusted), "\n")

cat("\nC2 - Significant pathways:\n")
cat("Raw p-value <= 0.05: ", nrow(c2_raw), "\n")
cat("Adjusted p-value <= 0.05: ", nrow(c2_adjusted), "\n")

# View top pathways for comparison
cat("\nTop Reactome Pathways:\n")
print(head(reactome_adjusted[, c("Description", "BgRatio", "GeneRatio", "pvalue", "p.adjust")]))

cat("\nTop C2:CP Pathways:\n")
print(head(c2cp_adjusted[, c("Description", "BgRatio", "GeneRatio", "pvalue", "p.adjust")]))

cat("\nTop C2 Pathways:\n")
print(head(c2_adjusted[, c("Description", "BgRatio", "GeneRatio", "pvalue", "p.adjust")]))



# Compare the Pathways Tested:
reactome_tested <- unique(reactome_df$Description)
c2cp_tested <- unique(c2cp_df$Description)
setdiff(reactome_tested, c2cp_tested) # Reactome pathways not in C2:CP
setdiff(c2cp_tested, reactome_tested) # C2:CP pathways not in Reactome

# Compare Gene Overlap:
reactome_genes <- reactome_df %>% filter(Description == "Signaling by the B Cell Receptor (BCR)") %>% pull(geneID)
c2cp_genes <- c2cp_df %>% filter(Description == "SA_B_CELL_RECEPTOR_COMPLEXES") %>% pull(geneID)
intersect(reactome_genes, c2cp_genes) # Shared genes

# Check Pathway Sizes:
reactome_sizes <- reactome_df %>% select(Description, BgRatio)
c2cp_sizes <- c2cp_df %>% select(Description, BgRatio)
reactome_sizes
c2cp_sizes


# checking Bgckground genes dynamic changed
# Compare Background Genes Used
# Check background genes indirectly from the BgRatio column
reactome_bg_genes <- reactome_df %>%
  mutate(BgCount = as.integer(sub("/.*", "", BgRatio))) %>%
  pull(BgCount)

c2cp_bg_genes <- c2cp_df %>%
  mutate(BgCount = as.integer(sub("/.*", "", BgRatio))) %>%
  pull(BgCount)

# Compare the number of background genes inferred
reactome_bg_count <- unique(reactome_bg_genes)
c2cp_bg_count <- unique(c2cp_bg_genes)

cat("Reactome inferred background gene count:", length(reactome_bg_count), "\n")
cat("C2:CP inferred background gene count:", length(c2cp_bg_count), "\n")


# checking input genes dynamic changed
# Check GeneRatio details for Reactome and C2:CP
reactome_gene_ratio <- reactome_df %>%
  mutate(
    GeneRatio_Numerator = as.integer(sub("/.*", "", GeneRatio)),
    GeneRatio_Denominator = as.integer(sub(".*/", "", GeneRatio))
  ) %>%
  select(Description, GeneRatio, GeneRatio_Numerator, GeneRatio_Denominator)

c2cp_gene_ratio <- c2cp_df %>%
  mutate(
    GeneRatio_Numerator = as.integer(sub("/.*", "", GeneRatio)),
    GeneRatio_Denominator = as.integer(sub(".*/", "", GeneRatio))
  ) %>%
  select(Description, GeneRatio, GeneRatio_Numerator, GeneRatio_Denominator)

# Display the results for Reactome
cat("Reactome GeneRatio details:\n")
print(head(reactome_gene_ratio))

# Display the results for C2:CP
cat("\nC2:CP GeneRatio details:\n")
print(head(c2cp_gene_ratio))

# Compare the number of unique GeneRatio denominators inferred
reactome_gene_count <- unique(reactome_gene_ratio$GeneRatio_Denominator)
c2cp_gene_count <- unique(c2cp_gene_ratio$GeneRatio_Denominator)

cat("\nReactome inferred unique GeneRatio denominators:", length(reactome_gene_count), "\n")
cat("C2:CP inferred unique GeneRatio denominators:", length(c2cp_gene_count), "\n")
