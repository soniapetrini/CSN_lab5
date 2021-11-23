
# PACKAGES

rm(list=ls())

library(igraph)
library(igraphdata)
library(ggplot2)
library(tidyr)
library(xtable)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/FIB/csn/lab/5/CSN_lab5")

# get nets
data("foodwebs")
data("macaque")
foodweb_net <- as.undirected(foodwebs[[9]], mode = "collapse")
macaque_net <- as.undirected(macaque, mode = "collapse")
karate_net <- graph.famous("Zachary")



## customized graph
merge_communities <- function(graphs, difficulty = length(graphs)+8 ) {
  
  # difficulty is an integer equal to the number of edges to add between communities
  if (difficulty < length(graphs)+8) {
    stop("Difficulty should be higher, or you could get an unconnected graph")
  }
  
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

# 2 version
merge_communities <- function(graphs, difficulty = NULL) {
  n_nodes_cum <- c(0)
  nodes_labs <- c(0)
  nodes <- list()
  EL_union <- c()
  t <- 1
  for (i in 1:length(graphs)) {
    g <- make_graph(graphs[i])
    n_nodes_cum[i+1] <- n_nodes_cum[i] + vcount(g) 
    nodes_labs[i+1] <- n_nodes_cum[i+1]-t
    # rename vertices
    labels <- seq(nodes_labs[i],nodes_labs[i+1])+1
    t %+=% 1
    vertex_attr(g) <- list(name = labels)
    nodes[[i]] <- V(g)
    EL  = get.edgelist(g)
    EL_union <- rbind(EL_union,EL)
  }
  GU = graph_from_edgelist(EL_union, directed=FALSE)
  
  # add some between-communities edges
  for (d in 1:difficulty) {
    r <- seq(1,length(nodes))
    index <- sample(r,2,replace = F)
    edge <- c(sample(nodes[[index[1]]],1),sample(nodes[[index[2]]],1))
    new_edge <- names(edge)
    print(new_edge)
    GU <- GU + edges(new_edge)
  }
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



## measures
get_measures <- function(graph, communities) {
  n <- vcount(graph)
  memberships <- membership(communities)
  clusters <- unique(memberships)
  
  # TPR
  nodes_in_tr_c <- c()
  i <- 1
  for (c in clusters) {
    members <- memberships[memberships==c]
    nodes <- names(members)
    subgraph <- induced_subgraph(graph,nodes,impl = "auto")
    nodes_in_triads <- unique(triangles(subgraph))
    nodes_in_tr_c[i] <- length(nodes_in_triads)   # number of nodes in triads
    i <- i + 1
  }
  TPR <- sum(nodes_in_tr_c)/n
  
  
  # Expansion
  is_edges_out = igraph::crossing(communities,graph)  # is the edge pointing out the cluster?
  fc_tot = 2*sum(is_edges_out)
  EXPANSION <- fc_tot/n
  
  
  # Conductance
  fc <- rep(0, length(clusters))
  mc <- rep(0, length(clusters))
  for (edge_index in 1:length(E(graph))) {
    edge <- list(origin=NULL,end=NULL)
    edge[c("origin","end")] <- V(graph)[ inc(edge_index) ]$name
    
    if (memberships[edge$origin] != memberships[edge$end]) {
      fc[memberships[edge$origin]] %+=% 1
      fc[memberships[edge$end]] %+=% 1}
    
    else {mc[memberships[edge$origin]] %+=% 1}
  }
  
  conductances <- fc/(2*mc+fc)
  weights <- table(memberships)/vcount(graph)
  CONDUCTANCE <- sum(conductances*weights)
  
  
  # Modularity
  MODULARITY <- modularity(communities)
  
  return(list("tpr"=TPR,"expansion"=EXPANSION,"conductance"=CONDUCTANCE,"modularity"= MODULARITY))
}




find_communities <- function(graph) {
  
  graph <- as.undirected(graph, mode='collapse')
  graph <- simplify(graph)
  vertex_attr(graph) <- list(name=V(graph))
  
  
  # summary
  N <- vcount(graph)
  E <- ecount(graph)
  N_dens <- N/E
  
  
  # edge.betweeness
  communities <- edge.betweenness.community(graph)
  n_clus_1 <- length(communities$membership %>% unique())
  measures1 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures1[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # fastgreedy.community
  communities <- fastgreedy.community(graph)
  n_clus_2 <- length(communities$membership %>% unique())
  measures2 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures2[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # label propagation
  communities <-  label.propagation.community(graph)
  n_clus_3 <- length(communities$membership %>% unique())
  measures3 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures3[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # leading eigenvector
  Isolated = which(degree(foodweb_net)==0)
  if (length(Isolated) != 0) {
    graph = delete.vertices(graph, Isolated)
  }
  
  arpack_defaults$maxiter = 1000000000
  communities <- leading.eigenvector.community(graph, options = arpack_defaults)
  n_clus_4 <- length(communities$membership %>% unique())
  measures4 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures4[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # louvain method (multilevel)
  communities <- cluster_louvain(graph)
  n_clus_5 <- length(communities$membership %>% unique())
  measures5 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures5[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # optimal clustering
  
  #?# not working on mac, package required
  
  
  # spinglass
  communities <- cluster_spinglass(graph)
  n_clus_7 <- length(communities$membership %>% unique())
  measures7 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures7[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # walktrap
  communities <- walktrap.community(graph)
  n_clus_8 <- length(communities$membership %>% unique())
  measures8 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures8[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  # infomap
  communities <- cluster_infomap(graph)
  n_clus_9 <- length(communities$membership %>% unique())
  measures9 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures9[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  df <- data.frame("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  df <- rbind(df,measures1,measures2,measures3,measures4,measures5,measures7,measures8,measures9)
  df$algorithm <- algorithms
  df$clusters <- c(n_clus_1,n_clus_2,n_clus_3,n_clus_4,n_clus_5,n_clus_7,n_clus_8,n_clus_9)
  
  return(df)
}


# STORAGE
notable_graphs <- c("Coxeter","Folkman","Herschel","Cubical")

my_easy_net <- merge_communities(notable_graphs,difficulty = 0)
my_hard_net <- merge_communities(notable_graphs,difficulty = 30)
plot(my_easy_net)
plot(my_hard_net)



nets <- list("macaque_net"=macaque_net,"karate_net"=karate_net,
             "my_easy_net"=my_easy_net,"my_hard_net"=my_hard_net,"foodwebs_net"=foodweb_net)


#for (i in 1:length(nets)) {
#  net <- nets[[i]]
#  net_name <- names(nets[i])
#  print(net_name)
#  df <- find_communities(net)
#  print(df)
#  #write.csv(df, paste(net_name,"df.csv",sep = "_"))
#}


df <- find_communities(net)
write.csv(df, paste(net_name,"df.csv",sep = "_"))






