library(sand)
library(igraphdata)

data(Ecoli.data)
g=graph.adjacency(regDB.adj)

d.coli <- degree(g)
par(mfrow=c(1,2))
hist(d.coli, xlab="Degree", ylab="Frequency", main="Degree Distribution")
dd.coli <- degree.distribution(g)
#vector of degree frequencies
d <- 1:max(d.coli)-1#u need this coz degree.distribution give exclusively freqs
ind <- (dd.coli != 0)
#log-log version
plot(d[ind], dd.coli[ind], log="xy", xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
 main="Log-Log Degree Distribution")

#neighbors degree vs degree to show neighbors paradox<------------------------------
a.nn.deg.coli <- graph.knn(g,V(g))$knn
plot(d.coli, a.nn.deg.coli, log="xy",
 col="goldenrod", xlab=c("Log Vertex Degree"),
 ylab=c("Log Average Neighbor Degree"))

l <- layout.kamada.kawai(g)
par(mfrow=c(1,2))
plot(g, layout=l, main="Hubs", vertex.label="",
 vertex.size=10 * sqrt(hub.score(g)$vector))#lot of hubs
#hub scores are the principal eigenvector of A*t(A) 
plot(g, layout=l, main="Authorities", vertex.label="",
 vertex.size=10 * sqrt(authority.score(g)$vector))#few authorities
#principal eigenvector of t(A)*A 

library(network)
library(sna)
A <- get.adjacency(g, sparse=FALSE)
g_sna <- network::as.network.matrix(A)
plot_radial<-function(measure){
sna::gplot.target(g_sna, measure, 
 circ.lab = FALSE, circ.col="skyblue",usearrows = FALSE,
 edge.col="darkgray")}
par(mfrow=c(2,2))
plot_radial(degree(g_sna))#few nodes have higher degree
plot_radial(betweenness(g_sna))# few nodes have higer betwenness
plot_radial(closeness(g_sna)) #no node has higher closseness
plot_radial(evcent(g_sna))#more dispersion than the other measures

#more on differences between igraph and sna closeness
igCl=igraph::closeness(g)
#Warning message:
#In igraph::closeness(g) :
#  At centrality.c:2784 :closeness centrality is not well-defined for disconnected graphs
snaCl=sna::closeness(g_sna)
sum(igCl==snaCl)
par(mfrow=c(1,2))
gplot.target(g_sna, igCl,circ.lab = FALSE, circ.col="skyblue",usearrows = FALSE,edge.col="darkgray")
gplot.target(g_sna, snaCl,circ.lab = FALSE, circ.col="skyblue",usearrows = FALSE,edge.col="darkgray")
#though measures are different, plots look the same
#maybe this is owed to the loooow values
max(igCl)
#[1] 5.171433e-05
#eitherway closeness has no sense for this net coz it is
# the inverse of the average length of the shortest paths TO/FROM ALL the other vertices

######################EDGES########################
g1=line.graph(g)
par(mfrow=c(1,2))
hist(edge.betweenness(g))
hist(betweenness(g1))
#no se parecen
#######################COHESION#######################
table(sapply(maximal.cliques(g), length))
# 1  2  3  4 
#35 73 36 18 
#density[0,1] provides a measure of how close a subgraph is to being a clique
graph.density(g)
#[1] 0.009459924

#a k-core is a maximal subgraph for which all vertex degrees are at least k
cores <- graph.coreness(g)
table(cores)
# 0  1  2  3  4 
#35 47 34 17 20
sna::gplot.target(g, cores, circ.lab = FALSE,
 circ.col="skyblue", usearrows = FALSE,
 vertex.col=cores, edge.col="darkgray")
detach("package:sna")
detach("package:network")

dyad.census(g_sna)
#     Mut Asym  Null
#[1,]  11  198 11419
#default=mutual dyads/asymmetric
reciprocity(g, mode="default")
#[1] 0.1
#ratio=mutual dyads/total number of edges
reciprocity(g, mode="ratio")
#[1] 0.05263158

#vector of occurences of each motif in the graph, ordered by isomorphism classes
graph.motifs()
moti=graph.motifs(g)
which(!is.na(moti))
#[1]  3  5  6  7  8  9 10 11 12 13 14 15 16

comps <- decompose.graph(g)#actual subgraphs not like components()
table(sapply(comps, vcount))#~components(g)$csize

#focus on larger component
gc <- decompose.graph(g)[[2]]
c(average.path.length(gc),average.path.length(g))
#[1] 2.645740 2.639857
c(diameter(gc),diameter(g))
#[1] 7 7
c(transitivity(gc),transitivity(g))
#[1] 0.08434783 0.08434783
# k-vertex-connected if the removal of any subset of vertices of cardinality < k
#leaves a subgraph that is connected
#vertex connectivity is the largest integer such that G is k-vertex- connected.
vertex.connectivity(gc)
#[1] 0
edge.connectivity(gc)
#[1] 0
# A single vertex that disconnects the graph is called a cut vertex, or articulation point.
cut.vertices <- articulation.points(gc)
length(cut.vertices)/vcount(gc)
#[1] 0.1607143
#maximal (weakly or strongly) connected components of the giant component
scc <- clusters(gc, mode=c("strong"))
table(scc$csize)
# 1  2  3  5 
#92  6  1  1 

#assortativity=selective linking among vertices~correlation coefficients
assortativity.degree(g)
#[1] -0.2921195

######################PARTITIONING########################
library(ape)

data(karate)
hierarchical <- fastgreedy.community(karate)#only for undirected
hierarchical1 <- cluster_fast_greedy(karate)
table(membership(hierarchical),membership(hierarchical1))
#gives the same partiotion, should be use cluster_fast_greedy?

k.lap <- graph.laplacian(karate)#spectral clustering
eig.anal <- eigen(k.lap)
f.vec <- eig.anal$vectors[, 33]
#how do u know when to stop spectral partioning of clusters?
leading<-cluster_leading_eigen(karate)
table(leading,ifelse(f.vec>0,1,0))
#leading_eigen != spectral clustering
#leading_eigen uses the leading non-negative eigenvector of the  modularity matrix 

optimal=cluster_optimal(karate)
#  At optimal_modularity.c:85 : GLPK is not available, Unimplemented function call
betw=cluster_edge_betweenness(karate)#weighted edge betweenness community detection might not make sense
label=cluster_label_prop(karate)
spin=cluster_spinglass(karate)#optimizing the an energy function.
info=cluster_infomap(karate)#minimizes trajectory of a random walker 
lou=cluster_louvain(karate)#hierarchical multi-level optimization of modularity
walk=cluster_walktrap(karate)#short random walks tend to stay in the same community

l <- layout.fruchterman.reingold(karate)
par(mfrow=c(2,4))
plot(hierarchical, karate,layout=l,main="hierarchical")
plot(karate,layout=l,vertex.color=membership(lou),main="louvain")
plot(karate,layout=l,vertex.color=ifelse(f.vec>0,"red","cyan"),main="spectral")
plot(karate,layout=l,vertex.color=membership(leading),main="leading vector")
plot(karate,layout=l,vertex.color=membership(spin),main="spinglass")
plot(karate,layout=l,vertex.color=membership(label),main="label propagation")
plot(karate,layout=l,vertex.color=membership(info),main="infomap")
plot(karate,layout=l,vertex.color=membership(walk),main="walktrap")

#spectral clustering
L <- graph.laplacian(g)
#L is a matrix L = D − A, with adjacency matrix A &
# D a diagonal matrix with the degree sequence 
eigstuff <- eigen(L)
summary(eigstuff$values)
# Length   Class    Mode 
#    153 complex complex 
####not amenable to this clustering?==============
