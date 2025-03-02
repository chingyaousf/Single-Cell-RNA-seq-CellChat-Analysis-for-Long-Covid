# Load required libraries
library(patchwork)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(dplyr)
library(Seurat)
library(tidyr)
library(openxlsx)
library(ggtext)

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(msigdbr) # this package provides access to the MSigDB gene sets
library(ReactomePA)
library(ggnewscale)


# using (max.cells.per.ident = 1000)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)



# Define the cell types to analyze separately
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
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
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Filter out RPL, RPS, and mitochondrial genes
  cell.markers <- cell.markers %>%
  filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene)  & !grepl("^MR-", gene))
  

  # Filter markers for each group
  longcovid_markers <- cell.markers %>% 
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>% 
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
  
  control_entrez <- bitr(control_markers$gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, 
                    Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis using the background
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Check if pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  

  # Define a new variable for the modified cell type name
  #cell_type_02 <- paste0(cell_type, "_02")  # Add "_02" to distinguish it from the original


  # Save the pathway enrichment results
  saveRDS(pathway_enrichment_comp, 
          file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, ".rds"))
  
  # Create a dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome.png"),
         clusterProfiler::dotplot(pathway_enrichment_comp,
                 showCategory = 15, 
                 label_format = 50,
                 color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")"))  +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 7, face = "bold", angle = 0, hjust = 1.1),  # Bold x-axis text
                axis.text.y = element_text(size = 7),  # Bold y-axis text
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
           scale_color_distiller(palette = "YlGnBu") +
           guides(fill = guide_colourbar(barwidth = 2, barheight = 10)),
         width = 5, 
         height = 6, 
         scale = 1.3,
         bg = "white",
         dpi = 600)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}




# without using (max.cells.per.ident = 1000)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)



# Define the cell types to analyze separately
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
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
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
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
  filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene)  & !grepl("^MR-", gene))
  

  # Filter markers for each group
  longcovid_markers <- cell.markers %>% 
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>% 
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
  
  control_entrez <- bitr(control_markers$gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, 
                    Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis using the background
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Check if pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  

  # Define a new variable for the modified cell type name
  cell_type_02 <- paste0(cell_type, "_02")  # Add "_02" to distinguish it from the original


  # Save the pathway enrichment results
  saveRDS(pathway_enrichment_comp, 
          file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type_02, ".rds"))
  
  # Create a dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type_02, "_LCvsControl_Reactome.png"),
         clusterProfiler::dotplot(pathway_enrichment_comp,
                 showCategory = 15, 
                 label_format = 50,
                 color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")"))  +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 7, face = "bold", angle = 0, hjust = 1.1),  # Bold x-axis text
                axis.text.y = element_text(size = 7),  # Bold y-axis text
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
           scale_color_distiller(palette = "YlGnBu") +
           guides(fill = guide_colourbar(barwidth = 2, barheight = 10)),
         width = 5, 
         height = 6, 
         scale = 1.3,
         bg = "white",
         dpi = 600)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}



# without using (max.cells.per.ident = 1000), NK only
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)



# Define the cell types to analyze separately
cell_types_to_analyze <- c("NK")

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
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
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
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
  filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene)  & !grepl("^MR-", gene))
  

  # Filter markers for each group
  longcovid_markers <- cell.markers %>% 
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>% 
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
  
  control_entrez <- bitr(control_markers$gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, 
                    Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis using the background
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Check if pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  

  # Define a new variable for the modified cell type name
  cell_type_02 <- paste0(cell_type, "_02")  # Add "_02" to distinguish it from the original


  # Save the pathway enrichment results
  saveRDS(pathway_enrichment_comp, 
          file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type_02, ".rds"))
  
  # Create a dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type_02, "_LCvsControl_Reactome_02.png"),
         clusterProfiler::dotplot(pathway_enrichment_comp,
                 showCategory = 15, 
                 label_format = 50,
                 color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")"))  +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 7, face = "bold", angle = 0, hjust = 1.1),  # Bold x-axis text
                axis.text.y = element_text(size = 7),  # Bold y-axis text
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
           scale_color_distiller(palette = "YlGnBu") +
           guides(fill = guide_colourbar(barwidth = 2, barheight = 10)),
         width = 5, 
         height = 6, 
         scale = 1.3,
         bg = "white",
         dpi = 600)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}








# using (max.cells.per.ident = 1000), fixing # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)



# Define the cell types to analyze separately
cell_types_to_analyze <- c("B_naive", "B_exhausted", "B_immature", "B_malignant", 
                           "B_memory", "CD16_mono", "CD4.Naive", "CD4.CM", "CD4.EM", 
                           "CD4.IL22", "CD4.Tfh", "CD8.Naive", "CD8.EM", "CD8.TE", 
                           "CD14_mono", "DC_cells", "HSC", "Lymphocytes", "MAIT", 
                           "NKT", "NK", "Plasma", "Platelets", "RBC", "Treg", "gdT", "pDC")

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
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
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  cell.markers <- cell.markers %>%
  filter(!grepl("^RPL", gene) & 
         !grepl("^RPS", gene) & 
         !grepl("^MT-", gene) & 
         !grepl("^MRPL|^MRPS", gene))
  

  # Filter markers for each group
  longcovid_markers <- cell.markers %>% 
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>% 
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
  
  control_entrez <- bitr(control_markers$gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, 
                    Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis using the background
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Check if pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  

  # Save all pathways to a CSV file for further inspection
  if (!is.null(pathway_enrichment_comp) && nrow(as.data.frame(pathway_enrichment_comp)) > 0) {
    write.csv(as.data.frame(pathway_enrichment_comp), 
              file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_03.csv"), 
              row.names = FALSE)
  }


  # Save the pathway enrichment results
  saveRDS(pathway_enrichment_comp, 
          file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_03.rds"))
  
  # Create a dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_03.png"),
         clusterProfiler::dotplot(pathway_enrichment_comp,
                 showCategory = 15, 
                 label_format = 50,
                 color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")"))  +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 7, face = "bold", angle = 0, hjust = 1.1),  # Bold x-axis text
                axis.text.y = element_text(size = 7),  # Bold y-axis text
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
           scale_color_distiller(palette = "YlGnBu") +
           guides(fill = guide_colourbar(barwidth = 2, barheight = 10)),
         width = 5, 
         height = 6, 
         scale = 1.3,
         bg = "white",
         dpi = 600)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}




# using (max.cells.per.ident = 1000), fixing # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
# (only cell types CD4.IL22, CD4.Tfh, CD8.TE, gdT)
# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Inspect the Seurat object metadata
head(seurat_object@meta.data)
unique(seurat_object@meta.data$disease)
unique(seurat_object@meta.data$CT_Sub_combined_02)



# Define the cell types to analyze separately
cell_types_to_analyze <- c("CD4.IL22", "CD4.Tfh", "CD8.TE", "gdT")
#cell_types_to_analyze <- c("gdT")

# Loop over each cell type for separate analysis
for (cell_type in cell_types_to_analyze) {
  
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
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set the identity to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Filter out unwanted genes (e.g., ribosomal, mitochondrial genes)
  cell.markers <- cell.markers %>%
  filter(!grepl("^RPL", gene) & 
         !grepl("^RPS", gene) & 
         !grepl("^MT-", gene) & 
         !grepl("^MRPL|^MRPS", gene))
  

  # Filter markers for each group
  longcovid_markers <- cell.markers %>% 
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>% 
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
  
  control_entrez <- bitr(control_markers$gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, 
                    Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis using the background
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Check if pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  

  # Save all pathways to a CSV file for further inspection
  if (!is.null(pathway_enrichment_comp) && nrow(as.data.frame(pathway_enrichment_comp)) > 0) {
    write.csv(as.data.frame(pathway_enrichment_comp), 
              file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_04.csv"), 
              row.names = FALSE)
  }


  # Save the pathway enrichment results
  saveRDS(pathway_enrichment_comp, 
          file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_04.rds"))
  
  # Create a dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_04.png"),
         clusterProfiler::dotplot(pathway_enrichment_comp,
                 showCategory = 10, 
                 label_format = 50,
                 color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")"))  +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 7, face = "bold", angle = 0, hjust = 1.1),  # Bold x-axis text
                axis.text.y = element_text(size = 7),  # Bold y-axis text
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
           scale_color_distiller(palette = "YlGnBu") +
           guides(fill = guide_colourbar(barwidth = 2, barheight = 10)),
         width = 5, 
         height = 6, 
         scale = 1.3,
         bg = "white",
         dpi = 600)
  
  message(paste("Pathway enrichment completed for cell type:", cell_type))
}





## change pathway text to "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)"
# Extract the data frame from compareClusterResult
pathway_results <- pathway_enrichment_comp@compareClusterResult

# Modify the pathway name
pathway_results$Description <- ifelse(
  pathway_results$Description == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
  "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)",
  pathway_results$Description
)

# Assign the modified data back to compareClusterResult
pathway_enrichment_comp@compareClusterResult <- pathway_results

# Verify if the name was changed
#unique(pathway_enrichment_comp@compareClusterResult$Description)

# Modify pathway name for highlighting
  pathway_enrichment_comp@compareClusterResult$Description <- ifelse(
    pathway_enrichment_comp@compareClusterResult$Description == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
    "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)",
    pathway_enrichment_comp@compareClusterResult$Description
  )
  
  # Generate dot plot with highlight
  plot <- clusterProfiler::dotplot(pathway_enrichment_comp, 
                                   showCategory = 10, 
                                   label_format = 50,
                                   color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")")) +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.9),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 9, face = "bold", angle = 0, hjust = 0.5),  
                axis.text.y = element_text(size = 9),  
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
          scale_color_distiller(palette = "YlGnBu") +
          guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) 

print(plot)

# Save the updated plot
ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_05.png"),
       plot, width = 5, height = 6, scale = 1.3, bg = "white", dpi = 600)

message(paste("Pathway enrichment completed for cell type:", cell_type))






# ===================================================
# Pathway Enrichment & Dot Plot Generation
# Long Covid vs Control (Seurat & Reactome)
# ## change pathway text to "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)"
# ===================================================
# Performs differential expression analysis, 
# pathway enrichment, and dot plot visualization 
# for selected cell types in Long Covid vs Control.
# Using (max.cells.per.ident = 1000)
# This script is optimized for analyzing specific cell types 
# (e.g., CD4.IL22, CD4.Tfh, CD8.TE, gdT) with pathway enrichment.
# Outputs the top 10 enriched Reactome pathways.

# Load required libraries
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define cell types to analyze
cell_types_to_analyze <- c("CD4.IL22", "CD4.Tfh", "CD8.TE", "gdT")

# Loop over each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset the current cell type
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Skip if no background genes are available
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set identity class to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run differential expression analysis
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Remove ribosomal, mitochondrial genes
  cell.markers <- cell.markers %>%
    filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene) & !grepl("^MRPL|^MRPS", gene))
  
  # Filter markers for each group
  longcovid_markers <- cell.markers %>%
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>%
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>%
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>%
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  control_entrez <- bitr(control_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Skip if no pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  
  # Extract the data frame from compareClusterResult
  pathway_results <- pathway_enrichment_comp@compareClusterResult

  # Modify pathway name for emphasis
  pathway_results$Description <- ifelse(
    pathway_results$Description == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
    "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)",
    pathway_results$Description
  )
  
  # Assign modified data back to compareClusterResult
  pathway_enrichment_comp@compareClusterResult <- pathway_results
  
  # Save results
  #write.csv(as.data.frame(pathway_enrichment_comp), 
  #         file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_05.csv"), 
  #         row.names = FALSE)
  
  #saveRDS(pathway_enrichment_comp, 
  #       file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_05.rds"))
  

  # Generate dot plot
  plot <- clusterProfiler::dotplot(pathway_enrichment_comp, 
                                   showCategory = 10, 
                                   label_format = 50,
                                   color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")")) +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 9, face = "bold", angle = 0, hjust = 0.5),  
                axis.text.y = element_text(size = 9, face = "bold"),  
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
          scale_color_distiller(palette = "YlGnBu") +
          guides(fill = guide_colourbar(barwidth = 2, barheight = 10))

  # Save the dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_05.png"),
         plot, width = 5, height = 5, scale = 1.3, bg = "white", dpi = 600)

  message(paste("Pathway enrichment completed for cell type:", cell_type))
}




# ===================================================
# Pathway Enrichment & Dot Plot Generation
# Long Covid vs Control (Seurat & Reactome)
# ## "without" change pathway text to "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)"
# ===================================================
# Performs differential expression analysis, 
# pathway enrichment, and dot plot visualization 
# for selected cell types in Long Covid vs Control.
# Using (max.cells.per.ident = 1000)
# This script is optimized for analyzing specific cell types 
# (e.g., CD4.IL22, CD4.Tfh, CD8.TE, gdT) with pathway enrichment.
# Outputs the top 10 enriched Reactome pathways.

# Load required libraries
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

# Load the Seurat object
seurat_object <- readRDS(file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds")

# Define cell types to analyze
cell_types_to_analyze <- c("CD4.IL22", "CD4.Tfh", "CD8.TE", "gdT")

# Loop over each cell type
for (cell_type in cell_types_to_analyze) {
  
  # Subset the current cell type
  current_cell_type_cells <- subset(seurat_object, CT_Sub_combined_02 == cell_type)
  
  # Extract all expressed genes as background
  expressed_genes <- rownames(current_cell_type_cells@assays$RNA@data)
  
  # Convert expressed genes to Entrez IDs for background
  background_entrez <- bitr(expressed_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Skip if no background genes are available
  if (is.null(background_entrez) || nrow(background_entrez) == 0) {
    message(paste("No background genes available for cell type:", cell_type))
    next
  }
  
  # Set identity class to 'disease' for comparison
  Idents(current_cell_type_cells) <- current_cell_type_cells$disease
  
  # Run differential expression analysis
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Remove ribosomal, mitochondrial genes
  cell.markers <- cell.markers %>%
    filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene) & !grepl("^MRPL|^MRPS", gene))
  
  # Filter markers for each group
  longcovid_markers <- cell.markers %>%
    filter(cluster == "LongCovid" & avg_log2FC > 0 & p_val_adj < 0.05) %>%
    dplyr::select(gene, avg_log2FC)
  
  control_markers <- cell.markers %>%
    filter(cluster == "Control" & avg_log2FC > 0 & p_val_adj < 0.05) %>%
    dplyr::select(gene, avg_log2FC)
  
  # Convert gene symbols to Entrez IDs
  longcovid_entrez <- bitr(longcovid_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  control_entrez <- bitr(control_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Create named lists for each condition
  gene_list <- list(LongCovid = longcovid_entrez$ENTREZID, Control = control_entrez$ENTREZID)
  
  # Perform pathway enrichment analysis
  pathway_enrichment_comp <- compareCluster(geneCluster = gene_list, 
                                            fun = "enrichPathway", 
                                            organism = "human", 
                                            readable = TRUE, 
                                            pvalueCutoff = 0.05, 
                                            universe = background_entrez$ENTREZID, 
                                            minGSSize = 10, 
                                            maxGSSize = 200)
  
  # Skip if no pathways are enriched
  if (is.null(pathway_enrichment_comp) || nrow(as.data.frame(pathway_enrichment_comp)) == 0) {
    message(paste("No pathways enriched for cell type:", cell_type))
    next
  }
  
  # Save results
  #write.csv(as.data.frame(pathway_enrichment_comp), 
            #file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_05.csv"), 
            #row.names = FALSE)
  
  #saveRDS(pathway_enrichment_comp, 
          #file = paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/pathway_enrichment_reactome_", cell_type, "_05.rds"))

  # Generate dot plot
  plot <- clusterProfiler::dotplot(pathway_enrichment_comp, 
                                   showCategory = 10, 
                                   label_format = 50,
                                   color = "p.adjust") +
          labs(x = "", size = "Gene Ratio", fill = "P-value", 
               title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")")) +
          theme_minimal() +
          theme(panel.border = element_rect(colour = "gray70", fill = NA), 
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                plot.title.position = "plot",
                axis.text.x = element_text(size = 9, face = "bold", angle = 0, hjust = 0.5),  
                axis.text.y = element_text(size = 9, face = "bold"),  
                legend.position = "right",
                legend.justification = "right", 
                legend.box = "vertical",
                legend.margin = margin(r = 0)) +
          scale_color_distiller(palette = "YlGnBu") +
          guides(fill = guide_colourbar(barwidth = 2, barheight = 10))

  # Save the dot plot
  ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_06.png"),
         plot, width = 5, height = 5, scale = 1.3, bg = "white", dpi = 600)

  message(paste("Pathway enrichment completed for cell type:", cell_type))
}










## still not working well, text overlapped
library(stringr)  # Ensure text wrapping works

# Unicode bold text for ICAM pathway
bold_pathway <- "𝗜𝗺𝗺𝘂𝗻𝗼𝗿𝗲𝗴𝘂𝗹𝗮𝘁𝗼𝗿𝘆 𝗶𝗻𝘁𝗲𝗿𝗮𝗰𝘁𝗶𝗼𝗻𝘀 𝗯𝗲𝘁𝘄𝗲𝗲𝗻 𝗮 𝗟𝘆𝗺𝗽𝗵𝗼𝗶𝗱 𝗮𝗻𝗱 𝗮 𝗻𝗼𝗻-𝗟𝘆𝗺𝗽𝗵𝗼𝗶𝗱 𝗰𝗲𝗹𝗹 (𝗜𝗖𝗔𝗠 𝘀𝗶𝗴𝗻𝗮𝗹𝗶𝗻𝗴 𝗽𝗮𝘁𝗵𝘄𝗮𝘆)"

# Modify pathway name and wrap text
pathway_results$Description <- ifelse(
  pathway_results$Description == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell (ICAM signaling pathway)",
  bold_pathway,
  str_wrap(pathway_results$Description, width = 60)  # Wrap long text
)

# Assign back to compareClusterResult
pathway_enrichment_comp@compareClusterResult <- pathway_results

# Generate dot plot
plot <- clusterProfiler::dotplot(pathway_enrichment_comp, 
                                 showCategory = 10, 
                                 label_format = 50,
                                 color = "p.adjust") +
        labs(x = "", size = "Gene Ratio", fill = "P-value", 
             title = paste("Over-represented Reactome pathways \nLong Covid vs Control (", cell_type, ")")) +
        theme_minimal() +
        theme(panel.border = element_rect(colour = "gray70", fill = NA), 
              plot.title = element_text(size = 12, face = "bold", hjust = 0.9),
              plot.title.position = "plot",
              axis.text.x = element_text(size = 9, face = "bold", angle = 0, hjust = 0.5),  
              axis.text.y = element_text(size = 9),  # No color condition here
              legend.position = "right",
              legend.justification = "right", 
              legend.box = "vertical",
              legend.margin = margin(r = 0)) +
        scale_color_distiller(palette = "YlGnBu") +
        guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) 

print(plot)

# Save the updated plot
ggsave(paste0("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis_04_exclude_MT_RPL_RPS_genes/Dotplot_", cell_type, "_LCvsControl_Reactome_06.png"),
       plot, width = 5, height = 6, scale = 1.3, bg = "white", dpi = 600)

message(paste("Pathway enrichment completed for cell type:", cell_type))












# Run FindAllMarkers for the comparison of LongCovid vs Control
  cell.markers <- FindAllMarkers(current_cell_type_cells, 
                                 assay = "RNA", 
                                 logfc.threshold = 0.25, 
                                 # max.cells.per.ident = 1000, 
                                 only.pos = TRUE)
  
  # Filter out RPL, RPS, and mitochondrial genes
  cell.markers <- cell.markers %>%
  filter(!grepl("^RPL", gene) & !grepl("^RPS", gene) & !grepl("^MT-", gene)  & !grepl("^MR-", gene))


# Identify any genes matching the unwanted patterns
unwanted_genes <- cell.markers %>%
  filter(grepl("^RPL", gene) | grepl("^RPS", gene) | grepl("^MT-", gene) | grepl("^MR-", gene))

# Display the unwanted genes
if (nrow(unwanted_genes) > 0) {
  print("Unfiltered unwanted genes found:")
  print(unwanted_genes)
} else {
  print("No unwanted genes found after filtering.")
}



unwanted_background_genes <- background_entrez %>%
  filter(grepl("^RPL", SYMBOL) | grepl("^RPS", SYMBOL) | grepl("^MT-", SYMBOL) | grepl("^MR-", SYMBOL))

if (nrow(unwanted_background_genes) > 0) {
  print("Unwanted genes found in the background gene set:")
  print(unwanted_background_genes)
} else {
  print("No unwanted genes in the background gene set.")
}



unwanted_entrez <- bitr(unwanted_genes$gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
print(unwanted_entrez)
