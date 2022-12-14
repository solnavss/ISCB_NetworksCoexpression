# ISCB LA 2022

# RIABIO 2022 Workshop: Analysis of biological and biomedical transcriptomic data using Artificial Intelligence (AI) and Machine Learning methods (ML) methods with R

## Session 2: Coexpression networks
#### Maribel Hernandez Rosales, CINVESTAV, Mexico
#### Marisol Navarro Miranda, CINVESTAV, Mexico


# igraph is a collection of network analysis tools.

![alt text](https://github.com/solnavss/ISCB_NetworksCoexpression/blob/main/images/igraph_R.png)

Lets look at [igraph funtions/tools](https://igraph.org/r/html/latest/).


``` r
setwd("/Users/solouli/Desktop/")

library(igraph)         # Network analysis library

```

### We generate a random network of 15 undirected nodes. We set a seed for random numbers generation.

``` r
set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)

plot(g,
     vertex.color="darkorchid",  # nodes colors
     vertex.size=20,             # size of nodes
     edge.color="black")         # edges colors
```

# We are going to use igraph functions to obtain some metrics:

Definition of the metrics can be foun in [Network Science](http://networksciencebook.com/) from Albert-László Barabási  

# 1) Vertexes

``` r
V(g)
```

# 2) Edges

``` r
E(g)
```

# 3) Degree

``` r
degree(g)

```
# 4) Distance matrix, paths

``` r
distances(g)
```

# 5) Connected components

``` r
is_connected(g)
count_components(g)
components(g)
```

# 6) Shortest path

``` r
all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)
```

# 7) Diameter

``` r
diameter(g)
```

# 8) Density

``` r
edge_density(g)
```

# 9) Betweenness centrality

``` r
betweenness(g)
```

# Write our network as graphml

``` r
write.graph(g,"/Users/solouli/Desktop/yeast_protein_interaction.graphml", format="graphml")
```

# Save our network as png

``` r
png("random_network.png", width = 300*10, height = 300*8,
    res = 300, units = "px")

set.seed(25032022)

plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

dev.off()
```

# > Study case: intronless genes

We are going to work in teams of 4 with a subsaet of intronless genes.

![alt text](https://github.com/solnavss/ISCB_NetworksCoexpression/blob/main/images/intronless.png)

You can find the published paper here: [Deciphering the Tissue-Specific Regulatory Role of Intronless Genes Across Cancers](hhttps://link.springer.com/chapter/10.1007/978-3-031-06220-9_18).

Lets look at [igraph funtions/tools](https://igraph.org/r/html/latest/).

![alt text](https://github.com/solnavss/ISCB_NetworksCoexpression/blob/main/images/paper_intronless.png)

# References

* Rsources

https://igraph.org/r/html/latest/
http://networksciencebook.com/chapter/1

https://kateto.net/2016/05/network-datasets/

https://kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
https://github.com/elaragon/R-igraph-Network-Workshop/blob/master/NetSciX%202016%20Workshop.R

* Data Bases

https://snap.stanford.edu/data/

https://networkrepository.com/

https://icon.colorado.edu/#!/networks

* Create a dataset

https://www.genecards.org/

https://string-db.org/cgi/input?sessionId=bvevnhF1MzII&input_page_show_search=on
