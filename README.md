# ISCB_NetworksCoexpression

# RIABIO 2022 Workshop: Analysis of biological and biomedical transcriptomic data using Artificial Intelligence (AI) and Machine Learning methods (ML) methods with R

## Session 2: Coexpression networks
#### Maribel Hernandez Rosales, CINVESTAV, Mexico
#### Marisol Navarro Miranda, CINVESTAV, Mexico

# igraph: igraph is a collection of network analysis tools with the emphasis on efficiency, portability and ease of use. igraph is open source and free. 


![alt text](https://github.com/solouli/Genomica_Computacional/blob/master/practica01/ncbi_01.png)


``` r
setwd("/Users/solouli/Desktop/")

library(igraph)         # Network analysis library

```

### We generate a random network of 15 undirected nodes. We set a seed for random numbers generation.

set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)

plot(g,
     vertex.color="darkorchid",  # nodes colors
     vertex.size=20,             # size of nodes
     edge.color="black")         # edges colors

# We are going to use igraph functions to obtain some metrics:

# 1) Vertexes
V(g)

# 2) Edges
E(g)

# 3) Degree
degree(g)

# 4) Distance matrix, paths
distances(g)

# 5) Connected components
is_connected(g)
count_components(g)
components(g)

# 6) Shortest path
all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)

# 7) Diameter
diameter(g)

# 8) Density
edge_density(g)

# 9) Degree distribution
degree_distribution(g)

# 10) Betweenness centrality
betweenness(g)

# 11) Clustering coefficient
transitivity(g)


#----------------------------------------------------------------------

# Case study 


# Descargamos la lista de aristas para la red de interacciones entre proteinas.

id <- "1O5NgBBewLPHGpKDFZfNmNVu0BY-Y-dDh"
yeast <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id))
head(yeast)

# Convertimos esta lista de aristas en una red.

yeast %>%
  graph_from_data_frame -> g

# Primer vistazo a la red
g

# Nodos de la red.
V(g)

get.data.frame(g, what="vertices") %>% head()

# Aristas de la red.
E(g)

get.data.frame(g, what="edges") %>% head()

plot(g,
     vertex.size=3,
     vertex.color="blue",
     vertex.label=NA,     # No queremos ver las etiquetas de los nodos.
     edge.arrow.size=0)   # No queremos ver las puntas de las flechitas.

# Usamos la función simplify() de igraph para quitar conexiones múltiples
# y autoconexiones. Usamos la función as.unidirected() para quitar las
# direcciones de las aristas.

g <- igraph::simplify(g)
g <- as.undirected(g)
g

# Extraemos el núcleo (core) de la red formado por los nodos que tienen
# grado mayor a 3. 

core <- coreness(g,mode="all")        # Vemos a qué core pertenecen los nodos.
core <- core[core>3]                  # Extraemos los nodos del core 3.
g <- induced_subgraph(g, names(core)) # Subgrafo inducido por los nodos en el vector core.
g

V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
get.data.frame(g,what="vertices") %>% head() 

# Visualizamos de modo que el tamaño de cada nodo sea proporcional a 
# su grado.

plot(g,
     vertex.size=V(g)$degree,
     vertex.color="blue",
     vertex.label=NA)

# Ajustamos el tamaño de los nodos.
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color="blue",
     vertex.label=NA)

# Lo mismo, pero el tamaño en función de la intermediación.

plot(g,
     vertex.size=log(V(g)$betweenness+1),
     vertex.color="blue",
     vertex.label=NA)

# Vemos cuáles son los nodos de mayor grado y mayor intermediación.

degree(g) %>% sort(decreasing=TRUE) %>% head(n=15)

get.data.frame(g,what="vertices") %>%
  as_tibble() %>%
  arrange(-degree) %>%
  head(n=15) %>%
  pull(name)

get.data.frame(g,what="vertices") %>%
  as_tibble() %>%
  arrange(-betweenness) %>%
  head(n=15) %>%
  pull(name)

# Visualizamos de nuevo. Pero esta vez salvamos.

png("network.png", width = 300*10, height = 300*8,
    res = 300, units = "px")

set.seed(25032022)
plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

dev.off()

# Guardamos la red en formato graphml para su uso futuro.

write.graph(g,"/Users/solouli/Desktop/yeast_protein_interaction.graphml", format="graphml")