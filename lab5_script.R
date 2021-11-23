
# PACKAGES

library(igraph)
library(igraphdata)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/FIB/csn/lab/5/CSN_lab5")

# get nets
data("foodwebs")
data("macaque")
foodweb_net <- as.undirected(foodwebs[[9]], mode = "collapse")
macaque_net <- as.undirected(macaque, mode = "collapse")
karate_net <- graph.famous("Zachary")


## customized graph

merge_communities <- function(graphs, difficulty = length(graphs)+4 ) {
  
  # difficulty is an integer equal to the number of edges to add between communities
  
  labels <- seq(1,1000,50)
  EL_union <- c()
  for (i in 1:length(graphs)) {
    g <- make_graph(graphs[i])
    # rename vertices
    possible_vertex_names <- seq(labels[i],labels[i+1])
    g_labels <- as.character(sample(possible_vertex_names,vcount(g)))
    vertex_attr(g) <- list(name = g_labels)
    EL  = get.edgelist(g)
    EL_union <- rbind(EL_union,EL)
  }
  GU = graph_from_edgelist(EL_union, directed=FALSE)
  
  # add some between-communities edges
  x <- combn(V(GU),2)
  my_edges <- x[,sample(seq(dim(x)[2]),difficulty,replace=F)]
  GU <- GU + edges(my_edges)
  return(GU)
}

`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))




# define globals
algorithms <- c("edge.betweenness",
                "fastgreedy",
                "label.propagation",
                "leading.eigenvector",
                "multilevel",
                "spinglass",
                "walktrap",
                "infomap")




# MEASURES

## Triangle Participation Ratio
get_TPR <- function(graph, communities) {
  n <- vcount(graph)
  memberships <- membership(communities)
  nodes_in_tr_c <- c()
  i <- 1
  for (m in unique(memberships)) {
    members <- memberships[memberships==m]
    nodes <- names(members)
    subgraph <- induced_subgraph(graph,nodes,impl = "auto")
    nodes_in_triads <- unique(triangles(subgraph))
    nodes_in_tr_c[i] <- length(nodes_in_triads)   # number of nodes in triads
    i <- i + 1
  }
  tpr <- sum(nodes_in_tr_c)/n
  return(tpr)
}



# Expansion new
get_expansion <-  function(graph, communities){
  
  memberships = membership(communities)
  n <- vcount(graph)
  is_edges_out = crossing(communities,graph)  # is the edge pointing out the cluster?
  fc_tot = 2*sum(is_edges_out)
  expansion <- fc_tot/n
  
  return(expansion)
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



## Conductance  fc/(2mc+fc)

get_conductance = function(graph, communities){
  
  edge_out = crossing(graph = graph,communities = communities)
  
}





# Conductance new
get_conductance <- function(graph,communities) {
  clusters <- unique(memberships)
  memberships <- membership(communities)
  vertex_attr(karate) <- list(name=V(karate))
  fc <- rep(0, length(clusters))
  mc <- rep(0, length(clusters))
  for (edge_index in 1:length(E(karate))) {
    edge <- list(origin=NULL,end=NULL)
    edge[c("origin","end")] <- V(karate)[ inc(edge_index) ]$name
    
    if (memberships[edge$origin] != memberships[edge$end]) {
      fc[memberships[edge$origin]] %+=% 1
      fc[memberships[edge$end]] %+=% 1}
    
    else {mc[memberships[edge$origin]] %+=% 1}
  }
  
  conductances <- fc/(2*mc+fc)
  weights <- table(memberships)/vcount(karate)
  return(sum(conductances*weights))
}



## Modularity (already implemented in iGraph)





# IMPLEMENTATION

find_communities <- function(graph) {
  
  graph <- as.undirected(graph, mode='collapse')
  graph <- simplify(graph)
  
  
  # edge.betweeness
  communities <- edge.betweenness.community(graph)
  print("edge.betweeness")
  n_clus_1 <- length(communities$membership %>% unique())
  mod_1 <- modularity(communities)
  tpr_1 <- get_TPR(graph, communities)
  exp_1 <- get_expansion(graph, communities)
  cond_1 <- get_conductance(graph, communities)
  
  # fastgreedy.community
  communities <- fastgreedy.community(graph)
  print("fastgreedy.community")
  n_clus_2 <- length(communities$membership %>% unique())
  mod_2 <- modularity(communities)
  tpr_2 <- get_TPR(graph, communities)
  exp_2 <- get_expansion(graph,communities)
  cond_2 <- get_conductance(graph, communities)
  
  # label propagation
  communities <-  cluster_label_prop(graph)
  print("label prop")
  n_clus_3 <- length(communities$membership %>% unique())
  mod_3 <- modularity(communities)
  tpr_3 <- get_TPR(graph, communities)  
  exp_3 <- get_expansion(graph, communities)
  cond_3 <- get_conductance(graph, communities)
  
  # leading eigenvector
  Isolated = which(degree(foodweb_net)==0)
  if (length(Isolated) != 0) {
    graph = delete.vertices(graph, Isolated)
  }
  
  arpack_defaults$maxiter = 1000000000
  communities <- leading.eigenvector.community(graph, options = arpack_defaults)
  print("leading")
  n_clus_4 <- length(communities$membership %>% unique())
  mod_4 <- modularity(communities)
  tpr_4 <- get_TPR(graph, communities)
  exp_4 <- get_expansion(graph, communities)
  cond_4 <- get_conductance(graph, communities)

  # louvain method (multilevel)
  communities <- cluster_louvain(graph)
  print("louvain")
  n_clus_5 <- length(communities$membership %>% unique())
  mod_5 <- modularity(communities)
  tpr_5 <- get_TPR(graph, communities)
  exp_5 <- get_expansion(graph,communities)
  cond_5 <- get_conductance(graph, communities)
  
  # optimal clustering

  #?# not working on mac, package required
  
  
  # spinglass
  communities <- cluster_spinglass(graph)
  print("spinglass")
  n_clus_7 <- length(communities$membership %>% unique())
  mod_7 <- modularity(communities)
  tpr_7 <- get_TPR(graph, communities) 
  exp_7 <- get_expansion(graph, communities)
  cond_7 <- get_conductance(graph, communities)
  
  # walktrap
  communities <- walktrap.community(graph)
  print("walktrap")
  n_clus_8 <- length(communities$membership %>% unique())
  mod_8 <- modularity(communities)
  tpr_8 <- get_TPR(graph, communities)
  exp_8 <- get_expansion(graph, communities)
  cond_8 <- get_conductance(graph, communities)
  
  # infomap
  communities <- cluster_infomap(graph)
  print("infomap")
  n_clus_9 <- length(communities$membership %>% unique())
  mod_9 <- modularity(communities)
  tpr_9 <- get_TPR(graph, communities)
  exp_9 <- get_expansion(graph, communities)
  cond_9 <- get_conductance(graph, communities)
  
  clusters <- c(n_clus_1,n_clus_2,n_clus_3,n_clus_4,n_clus_5,n_clus_7,n_clus_8,n_clus_9)
  modularities <- c(mod_1,mod_2,mod_3,mod_4,mod_5,mod_7,mod_8,mod_9)
  TPRs <- c(tpr_1,tpr_2,tpr_3,tpr_4,tpr_5,tpr_7,tpr_8,tpr_9)
  expansions <- c(exp_1,exp_2,exp_3,exp_4,exp_5,exp_7,exp_8,exp_9)
  conductances <- c(cond_1,cond_2,cond_3,cond_4,cond_5,cond_7,cond_8,cond_9)
  
  return(data.frame("algorithm"=algorithms,"clusters"=clusters,"modularity" = modularities, 
                    "TPR" = TPRs,"expansion" = expansions,"conductance" = conductances))
}



# STORAGE
notable_graphs <- c("Coxeter","Folkman","Herschel","Icosahedral","Cubical")
my_net <- merge_communities(notable_graphs[1:4],difficulty = 15)

nets <- list("macaque_net"=macaque_net,
             "karate_net"=karate_net,"my_net"=my_net)

#for (i in 1:length(nets)) {
#  net <- nets[[i]]
#  df <- find_communities(net)
#  write.csv(df, paste(names(nets[i]),"df.csv",sep = "_"))
#  print(df)
#}
#
#
#karate <- read.csv("karate_net_df.csv")
#
#len <- read.csv("algs_results_karate.csv")

