
library(igraph)
library(igraphdata)

# get nets
data("foodwebs")
foodweb_net <- as.undirected(foodwebs[[9]], mode = "collapse")

data("macaque")
macaque_net <- as.undirected(macaque, mode = "collapse")


# define globals
nets <- c(foodweb_net,macaque_net)
algorithms <- c("edge.betweenness.community",
                "fastgreedy.community",
                "label.propagation.community",
                "leading.eigenvector.community",
                "multilevel.community",
                "spinglass.community",
                "walktrap.community",
                "infomap.community")



# MEASURES

## Triangle Participation Ratio
get_TPR <- function(graph, communities) {
  memberships <- membership(communities)
  tpr_clusters <- c()
  i <- 1
  for (m in unique(memberships)) {
    members <- memberships[memberships==m]
    nodes <- names(members)
    subgraph <- induced_subgraph(graph,nodes,impl = "auto")
    nodes_in_triads <- unique(triangles(subgraph))
    tpr_clusters[i] <- length(nodes_in_triads)/vcount(subgraph)
    i <- i + 1
  }
  tpr <- mean(tpr_clusters)
  return(tpr)
}




# IMPLEMENTATION

find_communities <- function(graph) {
  
  graph <- as.undirected(graph, mode='collapse')
  
  # edge.betweeness
  communities <- edge.betweenness.community(graph)
  mod_1 <- modularity(communities)
  tpr_1 <- get_TPR(graph, communities)
  
  # fastgreedy.community
  communities <- fastgreedy.community(graph)
  mod_2 <- modularity(communities)
  tpr_2 <- get_TPR(graph, communities)
  
  # label propagation
  communities <-  cluster_label_prop(graph)
  mod_3 <- modularity(communities)
  tpr_3 <- get_TPR(graph, communities)  
  
  # leading eigenvector
  Isolated = which(degree(foodweb_net)==0)
  if (length(Isolated) != 0) {
    graph = delete.vertices(graph, Isolated)
  }
  
  arpack_defaults$maxiter = 1000000000
  communities <- leading.eigenvector.community(graph, options = arpack_defaults)
  mod_4 <- modularity(communities)
  tpr_4 <- get_TPR(graph, communities)   
  
  # louvain method (multilevel)
  communities <- cluster_louvain(graph)
  mod_5 <- modularity(communities)
  tpr_5 <- get_TPR(graph, communities)  
  
  # optimal clustering

  #◊# not working on mac, package required
  
  
  # spinglass
  communities <- cluster_spinglass(graph)
  mod_7 <- modularity(communities)
  tpr_7 <- get_TPR(graph, communities)  
  
  # walktrap
  communities <- walktrap.community(graph)
  mod_8 <- modularity(communities)
  tpr_8 <- get_TPR(graph, communities)  
  
  # infomap
  communities <- cluster_infomap(graph)
  mod_9 <- modularity(communities)
  tpr_9 <- get_TPR(graph, communities)  
  
  
  modularities <- c(mod_1,mod_2,mod_3,mod_4,mod_5,mod_7,mod_8,mod_9)
  TPRs <- c(tpr_1,tpr_2,tpr_3,tpr_4,tpr_5,tpr_7,tpr_8,tpr_9)
  
  return(data.frame("algorithm"=algorithms,"modularities" = modularities, "TPRs" = TPRs))

  
}



# STORAGE

measures_food <- find_communities(foodweb_net)

#for (net in nets) {
#  measures <- find_communities(foodweb_net)
#  measures_df <- data.frame("algorithm"=algorithms,"modularities" = measures$modularities,"TPRs" = measures$TPRs)
#}

write.csv(measures_food, "measures_foodwebs.csv")







