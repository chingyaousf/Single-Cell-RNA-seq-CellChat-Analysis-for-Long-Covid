# Single-Cell RNA-seq CellChat Analysis for Long Covid

## Overview

This repository contains the analysis of **cell-cell communication in Long Covid** using the **CellChat** package. The study focuses on identifying key **signaling pathways** and evaluating **interaction strength across various cell types**.

## Objectives

-   Identify and visualize cell-cell communication networks in **Long Covid vs. Control** samples.

-   Compare **signaling pathways** enriched in Long Covid.

-   Determine **key sender, receiver, mediator, and influencer** cell populations.

-   Highlight **inhibitory vs. activating interactions** across different immune cell types.

## Methodology

-   Single-cell RNA-seq data preprocessed using **Seurat**.

-   Cell-cell communication inferred using **CellChat**.

-   Pathway enrichment analysis using **Reactome**.

-   Visualization of interactions using **network graphs, river plots, and heatmaps**.

-   Statistical comparison of interaction strengths between **Long Covid and Control groups**.

## Dependencies

-   **R 4.3.0** or later

-   **Seurat (v5.1.0 / v4.3.0.1)**

-   **CellChat**

-   **ggplot2, dplyr, tidyr, ComplexHeatmap** (for visualization)

-   **ReactomePA, clusterProfiler** (for pathway enrichment)

## Data

-   Processed **Seurat object** containing **Long Covid and Control** samples.

-   CellChat objects stored in:

    -   `/home/cyang40/chingyao/long_covid_project/merged_samples_02/long_covid_control_cellchat/`

## Workflow

### 1. Data Preprocessing

-   **Filter & normalize** scRNA-seq data using Seurat.

-   **Cluster annotation** for major immune cell types.

### 2. CellChat Analysis

-   **Identify cell-cell interactions** based on ligand-receptor pairs.

-   **Compute interaction strength** and visualize networks.

-   **Compare Long Covid vs. Control** using differential interaction analysis.

### 3. Signaling Pathway Analysis

-   Perform **Reactome enrichment** on detected signaling pathways.

-   Highlight upregulated and downregulated pathways in Long Covid.

### 4. Visualization

-   **Network diagrams** showing cell-cell interactions.

-   **Bubble plots** highlighting key ligand-receptor pairs.

-   **Chord diagrams** depicting communication between cell types.

-   **Bar plots** summarizing interaction strength across conditions.

-   **Dot plots & heatmaps** summarizing pathway activation.

-   **Scatter plots** comparing signaling differences between Long Covid and Control.

-   **Violin plots** displaying gene expression of key signaling molecules.

-   **River plots** illustrating pathway communication.

```{=html}
<!-- -->
```
-   **Dot plots & heatmaps** summarizing pathway activation.

## Results

-   **Most connected cell types**: Monocytes, CD8 T cells, and DC cells.

-   **Key pathways** in Long Covid:

    -   **ICAM signaling** (ICAM2, ITGAL, ITGB2)

    -   **TGFβ pathway**

    -   **IFN-II pathway**

    -   **Resistin & IGF signaling**

-   **Comparison of inhibitory vs. activating interactions**.

    ![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/compare_interactions_27celltypes.png?raw=true){width="381"}

    ![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/distance_functional_similarity_LongCovid_Control_27celltypes_top15_fixed.png?raw=true){width="447"}

    ![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/informationFlow_selectedPathways.png?raw=true){width="701"}

    ![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/heatmap_incoming_comparison_27celltypes_02.png?raw=true)

    ![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/heatmap_outgoing_comparison_27celltypes_02.png?raw=true)

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/cellchat/ICAM_chord_CD4.Tfh_CD4.IL22_CD8.TE_gdT_comparison_02.png?raw=true){width="660"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Dotplot_CD4.IL22_LCvsControl_Reactome_03.png?raw=true){width="418"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Dotplot_CD4.Tfh_LCvsControl_Reactome_03.png?raw=true){width="400"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Dotplot_CD8.TE_LCvsControl_Reactome_03.png?raw=true){width="474"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Dotplot_gdT_LCvsControl_Reactome_03.png?raw=true){width="375"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Barplot_CD4.IL22_Reactome_FindMarkers.png?raw=true){width="323"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Barplot_CD4.Tfh_Reactome_FindMarkers.png?raw=true){width="372"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Barplot_CD8.TE_Reactome_FindMarkers.png?raw=true){width="375"}

![](https://github.com/chingyaousf/Single-Cell-RNA-seq-CellChat-Analysis-for-Long-Covid/blob/main/pathways/Barplot_gdT_Reactome_FindMarkers.png?raw=true){width="330"}
