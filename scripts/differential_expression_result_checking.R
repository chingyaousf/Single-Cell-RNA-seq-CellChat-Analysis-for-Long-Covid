# Load the pathway enrichment results (Reactome for NKT)
 
FindAllMarkers <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/Reactome_analysis/DEG_27celltypes/FindAllMarkers_CD8.TE.rds")

FindAllMarkers %>% arrange(desc(avg_log2FC)) %>% filter(avg_log2FC > 1,p_val_adj < 0.05)
