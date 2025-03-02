# Load necessary libraries
library(ggplot2)
library(igraph)
library(ggraph)

# Define the nodes (cell types)
cell_types <- c("Platelets", "Monocytes", "DC cells", "Lymphocytes", "Plasma cells", "Treg cells", "γδT cells", "CD4 T cells")

# Define the edges (pathways)
edges <- data.frame(
  from = c("Platelets", "Monocytes", "DC cells", "Plasma cells", "Monocytes", "Monocytes", 
           "Lymphocytes", "Plasma cells", "Monocytes", "DC cells"),
  to = c("Monocytes", "Monocytes", "DC cells", "γδT cells", "Lymphocytes", "DC cells",
         "Monocytes", "Treg cells", "Treg cells", "CD4 T cells"),
  pathway_type = c("Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", 
                   "Activating", "Activating", "Activating", "Activating")
)

# Create a graph object
g <- graph_from_data_frame(edges, directed = TRUE, vertices = cell_types)

# Plot the graph using ggraph
ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = pathway_type), arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'), show.legend = TRUE) +
  geom_node_point(size = 5, color = "skyblue") +
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 5) +
  scale_edge_color_manual(values = c("Inhibitory" = "red", "Activating" = "blue")) +
  theme_void() +
  ggtitle("Inhibitory and Activating Pathways between Cell Types")



# Define the nodes (cell types)
cell_types <- c("Platelets", "Monocytes", "DC cells", "Lymphocytes", "Plasma cells", "Treg cells", "γδT cells", "CD4 T cells")

# Define the edges (pathways)
edges <- data.frame(
  from = c("Platelets", "Monocytes", "DC cells", "Plasma cells", "Monocytes", "Monocytes", 
           "Lymphocytes", "Plasma cells", "Monocytes", "DC cells"),
  to = c("Monocytes", "Monocytes", "DC cells", "γδT cells", "Lymphocytes", "DC cells",
         "Monocytes", "Treg cells", "Treg cells", "CD4 T cells"),
  pathway_type = c("Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", 
                   "Activating", "Activating", "Activating", "Activating")
)

# Create a graph object
g <- graph_from_data_frame(edges, directed = TRUE, vertices = cell_types)

# Plot the graph using ggraph with a better layout
ggraph(g, layout = "kk") +  # 'kk' layout for better spread of nodes
  geom_edge_link(aes(color = pathway_type), arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'), show.legend = TRUE) +
  geom_node_point(size = 5, color = "skyblue") +
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 5) +
  scale_edge_color_manual(values = c("Inhibitory" = "red", "Activating" = "blue")) +
  theme_void() +
  ggtitle("Inhibitory and Activating Pathways between Cell Types")


# Define the nodes (cell types)
cell_types <- c("Platelets", "Monocytes", "DC cells", "Lymphocytes", "Plasma cells", "Treg cells", "γδT cells", "CD4 T cells", "gdT cells")

# Define the edges (pathways with multiple connections)
edges <- data.frame(
  from = c("Platelets", "Platelets", "Monocytes", "Monocytes", "Monocytes", 
           "DC cells", "DC cells", "Plasma cells", "Plasma cells", "Plasma cells", 
           "Lymphocytes", "Lymphocytes", "Monocytes", "Monocytes", "DC cells", 
           "Monocytes", "DC cells", "Monocytes", "DC cells"),
  to = c("Monocytes", "DC cells", "Monocytes", "DC cells", "Lymphocytes", 
         "Monocytes", "Lymphocytes", "γδT cells", "Monocytes", "DC cells", 
         "Monocytes", "DC cells", "Treg cells", "CD4 T cells", "Treg cells", 
         "gdT cells", "gdT cells", "Treg cells", "CD4 T cells"),
  pathway_type = c("Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", 
                   "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", "Inhibitory", 
                   "Activating", "Activating", "Activating", "Activating", "Activating", 
                   "Activating", "Activating", "Activating", "Activating")
)

# Create a graph object
g <- graph_from_data_frame(edges, directed = TRUE, vertices = cell_types)

# Plot the network using ggraph with clearer clustering
ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = pathway_type), arrow = arrow(length = unit(3, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(size = 5, color = "skyblue") +
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 5) +
  scale_edge_color_manual(values = c("Inhibitory" = "red", "Activating" = "blue")) +
  theme_void() +
  ggtitle("Inhibitory and Activating Pathways between Cell Types")
