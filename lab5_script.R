
# PACKAGES

library(igraph)
library(igraphdata)
library(ggplot2)
library(tidyr)
library(xtable)

defaultW <- getOption("warn") 
options(warn = -1) 

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/FIB/csn/lab/5/CSN_lab5")

# get nets
data("foodwebs")
data("macaque")
data("rfid")
foodweb_net <- as.undirected(foodwebs[[9]], mode = "collapse")
macaque_net <- as.undirected(macaque, mode = "collapse")
rfid_net <- as.undirected(rfid, mode = "collapse")
karate_net <- graph.famous("Zachary")



## customized graph (communities with ER graphs)
merge_communities_ER <- function(comm_sizes,ERp, difficulty = NULL) {
  n_communities <- length(comm_sizes)
  vcount_cum <- c(0)
  nodes <- list()
  EL_union <- c()
  t <- 1
  for (i in 1:n_communities) {
    g <- erdos.renyi.game(comm_sizes[i],ERp)
    vcount_cum[i+1] <- vcount_cum[i] + vcount(g) 
    labels <- seq(vcount_cum[i]+ 1,vcount_cum[i+1])  
    vertex_attr(g) <- list(name = labels)
    nodes[[i]] <- V(g)
    EL_union <- rbind(EL_union,get.edgelist(g))
  }
  GU = graph_from_edgelist(EL_union, directed=FALSE)
  
  # add some between-communities edges
  difficulty = difficulty*5
  for (d in 1:difficulty) {
    clusters <- sample(seq(1,n_communities),2,replace = F)
    edge <- c(sample(nodes[[clusters[1]]],1),sample(nodes[[clusters[2]]],1))
    new_edge <- names(edge)
    GU <- GU + edges(new_edge)
  }
  plot(GU)
  return(GU)
}

# customized graph 3 version (communities with complete graphs)
merge_communities_full <- function(comm_sizes, difficulty = NULL) {
  n_communities <- length(comm_sizes)
  vcount_cum <- c(0)
  nodes <- list()
  EL_union <- c()
  t <- 1
  for (i in 1:n_communities) {
    g <- graph.full(comm_sizes[i])
    vcount_cum[i+1] <- vcount_cum[i] + vcount(g) 
    labels <- seq(vcount_cum[i]+ 1,vcount_cum[i+1])    # rename vertices with sequential numbers 
    vertex_attr(g) <- list(name = labels)              # to have minimum bridges
    nodes[[i]] <- V(g)
    EL  = get.edgelist(g)
    EL_union <- rbind(EL_union,EL) }
  GU = graph_from_edgelist(EL_union, directed=FALSE)   # connect the communities
  
  difficulty = difficulty*5
  for (d in 1:difficulty) {
    clusters <- sample(seq(1,n_communities),2,replace = F)
    edge <- c(sample(nodes[[clusters[1]]],1),sample(nodes[[clusters[2]]],1))
    new_edge <- names(edge)
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


plot(karate_net)

get_stats <- function(graph) {
  # summary
  N <- vcount(graph)
  E <- ecount(graph)
  E_dens <- E/N
  k_mean <- mean(degree(graph))
  graph_stats <- c(N,E,E_dens,k_mean)
  return(graph_stats)
}



find_communities <- function(graph) {
  
  graph <- as.undirected(graph, mode='collapse')
  graph <- simplify(graph)
  vertex_attr(graph) <- list(name=V(graph))
  communities <-  label.propagation.community(graph)
  
  
  ptm <- Sys.time()
  # edge.betweeness
  communities <- edge.betweenness.community(graph)
  print("edge.betweeness")
  time_elapsed1 <- Sys.time() - ptm
  print(time_elapsed1)
  n_clus_1 <- length(communities$membership %>% unique())
  measures1 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures1[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # fastgreedy.community
  communities <- fastgreedy.community(graph)
  print("fastgreedy.community")
  time_elapsed2 <- Sys.time() - ptm
  print(time_elapsed2)
  n_clus_2 <- length(communities$membership %>% unique())
  measures2 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures2[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # label propagation
  communities <-  label.propagation.community(graph)
  print("label")
  time_elapsed3 <- Sys.time() - ptm
  print(time_elapsed3)
  n_clus_3 <- length(communities$membership %>% unique())
  measures3 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures3[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()          # timing leading.eigenvector.community()
  
  # leading eigenvector
  print("eigenvector")
  Isolated = which(degree(graph)==0)
  if (length(Isolated) != 0) {
    graph = delete.vertices(graph, Isolated)
  }
  
  arpack_defaults$maxiter = 1000000000
  communities <- leading.eigenvector.community(graph, options = arpack_defaults)
  time_elapsed4 <- Sys.time() - ptm
  print(time_elapsed4)
  n_clus_4 <- length(communities$membership %>% unique())
  measures4 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures4[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # louvain method (multilevel)
  communities <- cluster_louvain(graph)
  print("louvain")
  time_elapsed5 <- Sys.time() - ptm
  print(time_elapsed5)
  n_clus_5 <- length(communities$membership %>% unique())
  measures5 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures5[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # spinglass
  communities <- cluster_spinglass(graph)
  print("spinglass")
  time_elapsed7 <- Sys.time() - ptm
  print(time_elapsed7)
  n_clus_7 <- length(communities$membership %>% unique())
  measures7 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures7[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # walktrap
  communities <- walktrap.community(graph)
  print("walktrap")
  time_elapsed8 <- Sys.time() - ptm
  print(time_elapsed8)
  n_clus_8 <- length(communities$membership %>% unique())
  measures8 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures8[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  ptm <- Sys.time()
  # infomap
  communities <- cluster_infomap(graph)
  print("infomap")
  time_elapsed9 <- Sys.time() - ptm
  print(time_elapsed9)
  n_clus_9 <- length(communities$membership %>% unique())
  
  measures9 <- list("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  measures9[c("tpr","expansion","conductance","modularity")] <- get_measures(graph,communities)
  
  df <- data.frame("tpr"=NULL,"expansion"=NULL,"conductance"=NULL,"modularity"= NULL)
  df <- rbind(df,measures1,measures2,measures3,measures4,measures5,measures7,measures8,measures9)
  df$algorithm <- algorithms
  df$clusters <- c(n_clus_1,n_clus_2,n_clus_3,n_clus_4,n_clus_5,n_clus_7,n_clus_8,n_clus_9)
  df$time <- c(time_elapsed1,time_elapsed2,time_elapsed3,time_elapsed4,time_elapsed5,
               time_elapsed7,time_elapsed8,time_elapsed9)
  return(df)
}


options(warn = defaultW)



### build from complete graphs  
comm_sizes <- c(94,53,25,62)
full_net <- merge_communities_full(comm_sizes,difficulty = 100)
#plot(full_net)
#save(full_net,file="full_net.rda")


## build from ER graphs  
comm_sizes <- c(94,53,25,62)
ERp <- 0.6
ER_net <- merge_communities_ER(comm_sizes,ERp,difficulty = 100)
plot(ER_net)

## get results

nets <- list("macaque_net"=macaque_net,"rfid_net"=rfid_net,
             "full_net"=full_net,"ER_net"=ER_net)


# summary table
sum_tab <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in 1:length(nets)){
  net_stats <- lapply(nets[i],FUN=get_stats)
  net_stats <- append(unlist(net_stats),names(nets[i]))
  sum_tab <- rbind(sum_tab,net_stats)
}
colnames(sum_tab) <- c("N","E","E_dens","k_mean","network")




for (i in 1:length(nets)) {
  net <- nets[[i]]
  net_name <- names(nets[i])
  print(net_name)
  df <- find_communities(net)
  #write.csv(df, paste(net_name,"df.csv",sep = "_"))
}



# TASK 2
## Wikipedia gml

### Loading data and Initial exploration

require(tm)
require(ggplot2)
require(wordcloud)
require(SnowballC)
require(tidytext)
require(RColorBrewer)

wiki = read.graph("wikipedia.gml", format="gml")
uwiki = as.undirected(wiki,mode = "collapse") # working with undirected version

vcount(wiki); ecount(wiki)
edge_density(uwiki,loops = F)

d = degree(uwiki,mode = "all")

as.data.frame(table(d)) %>%
  ggplot(aes(x = as.numeric(d),y = Freq)) +
  geom_point() +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  labs(title = "Summary Wikipedia dataset",subtitle = "Degree distribution and descriptives",
       x="Degree",y="Number of vertices") +
  ggplot2::annotate("label",label="Number of vertices: 27475\nNumber of edges: 85729\nEdge density: 0.0002",
                    x = 3,y=10)

### Community detection


cluster = function(graph, clust_alg, complexity = 5){
  
  c = clust_alg(graph)
  memberships = membership(c)
  
  singletones = which(table(c$membership) < complexity)
  
  keep = V(graph)[!(c$membership %in% singletones)]
  
  graph = induced.subgraph(graph,keep)
  
  c = clust_alg(graph) #clean graph
  
  n <- vcount(graph)
  
  is_edges_out = igraph::crossing(c,graph)  # is the edge pointing out the cluster?
  fc_tot = 2*sum(is_edges_out)
  
  EXPANSION <- fc_tot/n
  modularity = modularity(c)
  n_groups = length(unique(c$membership))
  
  return(data.frame(algorithm = c$algorithm, n_groups = n_groups, expansion = EXPANSION ,modularity = modularity, row.names = NULL))
}


wiki_results = rbind(cluster(uwiki, fastgreedy.community),
                     cluster(uwiki, walktrap.community),
                     cluster(uwiki, label.propagation.community),
                     cluster(uwiki, infomap.community))

## Best clust algorithm chosen: fast greedy

c = fastgreedy.community(uwiki)
memberships = membership(c)

singletones = which(table(c$membership) < 5)

keep = V(graph)[!(c$membership %in% singletones)]

graph = induced.subgraph(uwiki,keep)

c = fastgreedy.community(graph)


## Result visualization of the 

cloud_importance = function(group){
  
  meaningless = c("of","and","in","the","for", "a", "to", "The", "List")
  x = sapply(group, function(x) strsplit(x, split = " "),simplify = T,USE.NAMES = F)
  x = unlist(x)
  
  tab = as.data.frame(table(x))
  tab = tab[order(tab$Freq,decreasing = T),]
  tab = tab[!(tab$x %in% meaningless),]
  
  wordcloud(words = tab$x, tab$Freq, max.words = 30,random.order=FALSE,
            rot.per=0.35,colors=brewer.pal(8, "Dark2"))
}

# Plotting grid: cloud words 

par(mfrow=c(3,4))

for(i in 1:12){
  group = V(graph)[c$membership==i]$label
  cloud_importance(group)
}







