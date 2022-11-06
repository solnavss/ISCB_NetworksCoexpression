############################################################################
# ISCB LA 2022
# RIABIO 2022 Workshop: Analysis of biological and biomedical transcriptomic
# data using Artificial Intelligence (AI) and Machine Learning methods (ML) 
# methods with R
# Session 2: Coexpression networks
# Maribel Hernandez Rosales, CINVESTAV, Mexico
# Marisol Navarro Miranda, CINVESTAV, Mexico
############################################################################

#setwd('~/Desktop')
#install.packages(igraph)

library(igraph)

set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)
set.seed(25032022)

plot(g,
     vertex.color="darkorchid",  # color de los nodos
     vertex.size=20,             # tamaño de los nodos
     edge.color="black")         # color de las aristas

# 1) Vertix

V(g)

# 2) Edges

E(g)

# 3) El grado de los nodos | Degree

degree(g)

# 4) Matriz de distancias, caminos | Distance matrix, paths

distances(g)

# 5) Componentes conexos | Connected components

is_connected(g)
count_components(g)
components(g)

# 6) Camino mas corto | Shortest path

all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)

# 7) Diametro | Diameter

diameter(g)

# 8) Densidad | Density

edge_density(g)

# 9) Distribución de grado | Degree distribution

degree_distribution(g)

# 10) Centralidad,coeficiente de intermediación | Betweenness centrality

betweenness(g)

# 11) Coeficiente de clustering | Clustering coefficient

transitivity(g)

############################################################################
# Evolutionary Perspective and Expression Analysis of Intronless Genes 
# Highlight the Conservation of Their Regulatory Role
############################################################################

setwd("~/Desktop/")
library(igraph)

# HEALTHY

healthy <- read.table(file = '~/Desktop/coexp_h_5.tsv', sep = '\t', header = TRUE)
head(healthy)
class(healthy)

healthy_graph <- graph.data.frame(healthy, directed = FALSE)
head(healthy_graph)
class(healthy_graph)

plot(healthy_graph)

plot(healthy_graph,
     vertex.color="darkorchid",  # nodes color
     vertex.size=20,             # nodes size
     edge.color="black")         # edges color

# TUMOR

tumor <- read.table(file = '~/Desktop/coexp_t_5.tsv', sep = '\t', header = TRUE)
head(tumor)
class(tumor)

tumor_graph <- graph.data.frame(tumor, directed = FALSE)
head(tumor_graph)
class(tumor_graph)

plot(tumor_graph)

plot(tumor_graph,
     vertex.color="salmon",  # nodes color
     vertex.size=20,             # nodes size
     edge.color="black")         # edges color

############################################################################