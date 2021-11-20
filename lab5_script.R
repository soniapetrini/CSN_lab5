
# PACKAGES

library(igraph)
library(igraphdata)

rm(list=ls())

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
  w <- c()
  for (m in unique(memberships)) {
    members <- memberships[memberships==m]
    nodes <- names(members)
    subgraph <- induced_subgraph(graph,nodes,impl = "auto")
    sub_nodes <- vcount(subgraph)
    w[i] <- sub_nodes
    nodes_in_triads <- unique(triangles(subgraph))
    tpr_clusters[i] <- length(nodes_in_triads)/sub_nodes
    i <- i + 1
  }
  tpr <- weighted.mean(tpr_clusters, w)
  return(tpr)
}

## Expansion 2 (fc/nc)

get_expansion_2 = function(graph, communities){
  
  memberships = membership(communities)
  expansion = c()
  vertex_out = crossing(memberships) # is the edge pointing out the cluster?
  in_out = table(vertex_out)
  fc = in_out["TRUE"]
  nc = sum(in_out)
  
  return(fc/nc)
}

## Expansion
get_expansion <- function(graph,communities) {
  memberships <- membership(communities)
  list <- c()
  i <- 1
  for (edge_index in 1:length(E(graph))) {
    edge_nodes <- V(graph)[ inc(edge_index) ]$name
    if (memberships[edge_nodes[1]] != memberships[edge_nodes[2]]) {
      list[i] <- memberships[edge_nodes[1]]
      list[i+1] <- memberships[edge_nodes[2]]
      i <- i + 2
    }
  }
  expansions <- table(list)/table(memberships)
  expansion <- mean(expansions)
  return(expansion)
}

## Conductance

get_conductance = function(graph, communities){
  
  edge_out = crossing(graph = graph,communities = communities)
  
  
  
}

## Modularity (already implemented in iGraph)






# IMPLEMENTATION

find_communities <- function(graph) {
  
  graph <- as.undirected(graph, mode='collapse')
  graph <- simplify(graph)
  
  # edge.betweeness
  communities <- edge.betweenness.community(graph)
  mod_1 <- modularity(communities)
  tpr_1 <- get_TPR(graph, communities)
  exp_1 <- get_expansion(graph, communities)
  
  # fastgreedy.community
  communities <- fastgreedy.community(graph)
  mod_2 <- modularity(communities)
  tpr_2 <- get_TPR(graph, communities)
  exp_2 <- get_expansion(graph,communities)
  
  # label propagation
  communities <-  cluster_label_prop(graph)
  mod_3 <- modularity(communities)
  tpr_3 <- get_TPR(graph, communities)  
  exp_3 <- get_expansion(graph, communities)
  
  # leading eigenvector
  Isolated = which(degree(foodweb_net)==0)
  if (length(Isolated) != 0) {
    graph = delete.vertices(graph, Isolated)
  }
  
  arpack_defaults$maxiter = 1000000000
  communities <- leading.eigenvector.community(graph, options = arpack_defaults)
  mod_4 <- modularity(communities)
  tpr_4 <- get_TPR(graph, communities)
  exp_4 <- get_expansion(graph, communities)

  # louvain method (multilevel)
  communities <- cluster_louvain(graph)
  mod_5 <- modularity(communities)
  tpr_5 <- get_TPR(graph, communities)
  exp_5 <- get_expansion(graph,communities)
  
  # optimal clustering

  #â# not working on mac, package required
  
  
  # spinglass
  communities <- cluster_spinglass(graph)
  mod_7 <- modularity(communities)
  tpr_7 <- get_TPR(graph, communities) 
  exp_7 <- get_expansion(graph, communities)
  
  # walktrap
  communities <- walktrap.community(graph)
  mod_8 <- modularity(communities)
  tpr_8 <- get_TPR(graph, communities)
  exp_8 <- get_expansion(graph, communities)
  
  # infomap
  communities <- cluster_infomap(graph)
  mod_9 <- modularity(communities)
  tpr_9 <- get_TPR(graph, communities)
  exp_9 <- get_expansion(graph, communities)
  
  
  modularities <- c(mod_1,mod_2,mod_3,mod_4,mod_5,mod_7,mod_8,mod_9)
  TPRs <- c(tpr_1,tpr_2,tpr_3,tpr_4,tpr_5,tpr_7,tpr_8,tpr_9)
  expansions <- c(exp_1,exp_2,exp_3,exp_4,exp_5,exp_7,exp_8,exp_9)
  
  return(data.frame("algorithm"=algorithms,"modularity" = modularities, 
                    "TPR" = TPRs,"expansion" = expansions))
}



# STORAGE

measures_food <- find_communities(foodweb_net)

#for (net in nets) {
#  measures <- find_communities(foodweb_net)
#  measures_df <- data.frame("algorithm"=algorithms,"modularities" = measures$modularities,"TPRs" = measures$TPRs)
#}

write.csv(measures_food, "measures_foodwebs.csv")





