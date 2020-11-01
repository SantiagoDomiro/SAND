library(sand)
data(karate)

#######################DEGREE########################
hist(degree(karate), col="lightblue", xlim=c(0, 50),
 xlab="Vertex Degree", ylab="Frequency", main="") 

#strength=sum weights of edges incident to a vertex
#strength distribution=weighted degree distribution
distributionhist(graph.strength(karate), col="pink",
 xlab="Vertex Strength", ylab="Frequency", main="")

library(igraphdata)
data(yeast)
summary(yeast)
d.yeast <- degree(yeast)
#degree per gene
hist(d.yeast,col="blue", xlab="Degree", ylab="Frequency", main="Degree Distribution")

dd.yeast <- degree.distribution(yeast)
#vector of degree frequencies
d <- 1:max(d.yeast)-1
ind <- (dd.yeast != 0)
#log-log version
plot(d[ind], dd.yeast[ind], log="xy", col="blue",
 xlab=c("Log-Degree"), ylab=c("Log-Intensity"),
 main="Log-Log Degree Distribution")

a.nn.deg.yeast <- graph.knn(yeast,V(yeast))$knn
#average nearest neighbor degree per gene
plot(d.yeast, a.nn.deg.yeast, log="xy",
 col="goldenrod", xlab=c("Log Vertex Degree"),
 ylab=c("Log Average Neighbor Degree"))
#do nodes link with nodes of the same degree?

#radial layout to display vertex centralities
library(network)
#network class = range of relational data types, supports arbitrary attributes.
library(sna)
# Tools for Social Network Analysis
A <- get.adjacency(karate, sparse=FALSE)
g <- network::as.network.matrix(A)
plot_radial<-function(measure){
sna::gplot.target(g, measure, 
 circ.lab = FALSE, circ.col="skyblue",usearrows = FALSE,
 vertex.col=c("blue", rep("red", 32), "yellow"),
 edge.col="darkgray")}

######################OTHER CENTRALITY MEASURES########################
#hubs are not top for other centralities
par(mfrow=c(2,2))
#here ur using sna functions
plot_radial(degree(g))
plot_radial(betweenness(g))# vertices that sit on many paths is more critical communication
plot_radial(closeness(g)) #a vertex is central if it is ‘close’ to many other 
#def, usage and results are not equal with igraph function
plot_radial(evcent(g))
#evcent=eigen centrality= the more central the neighbors, the more central the vertex is
#result has no names, so u plot this not evcent(g)$vector as the book says

######################FOR DIRECTED NETS########################
l <- layout.kamada.kawai(aidsblog)
par(mfrow=c(1,2))
plot(aidsblog, layout=l, main="Hubs", vertex.label="",
 vertex.size=10 * sqrt(hub.score(aidsblog)$vector))
#hub scores are the principal eigenvector of A*t(A) 
plot(aidsblog, layout=l, main="Authorities", vertex.label="",
 vertex.size=10 * sqrt(authority.score(aidsblog)$vector))
#authorities are the principal eigenvector of t(A)*A
#they are not the same

######################EDGE CENTRALITY########################
eb <- edge.betweenness(karate)
E(karate)[order(eb, decreasing=T)[1:3]]

# to apply vertex centrality measures change vertices to edges, and edges to vertices
plot(karate)
plot(line.graph(karate))

######################COHESION########################
#A census of cliques can provide some sense of how structured a graph is.
table(sapply(cliques(karate), length))
#get the largest cliques
cliques(karate)[sapply(cliques(karate), length) == 5]
#loose mamushka cliques
table(sapply(maximal.cliques(karate), length))

#a k-core is a maximal subgraph for which all vertex degrees are at least k
cores <- graph.coreness(karate)
#core per node
sna::gplot.target(g, cores, circ.lab = FALSE,
 circ.col="skyblue", usearrows = FALSE,
 vertex.col=cores, edge.col="darkgray")
detach("package:network")
detach("package:sna")

#A census of dyad state can yield insight into the connectivity in the graph
aidsblog <- simplify(aidsblog)
dyad.census(aidsblog)
#$mut = mutual (two directed edges).
# [1] 3
# $asym = asymmetric (one directed edge)
# [1] 177
# $null = no directed edges
#[1] 10405

#ego net for the instructor and the admin in karate
ego.instr <- induced.subgraph(karate,neighborhood(karate, 1, 1)[[1]])
ego.admin <- induced.subgraph(karate, neighborhood(karate, 1, 34)[[1]])
#density[0,1] provides a measure of how close a subgraph is to being a clique
graph.density(karate)
#[1] 0.1390374
graph.density(ego.instr)
# [1] 0.25
graph.density(ego.admin)
#[1] 0.2091503

#transitivity=clustering
transitivity(karate)
#[1] 0.2556818
transitivity(karate, "local", vids=c(1,34))
#[1] 0.1500000 0.1102941

#for directed graphs  reciprocity has to definition: default=mutual dyads/asymmetric
#ratio=mutual dyads/total number of edges
reciprocity(aidsblog, mode="default")
# [1] 0.03278689
reciprocity(aidsblog, mode="ratio")
#[1] 0.01666667

#components
comps <- decompose.graph(yeast)#acutal subgraphs not like components()
table(sapply(comps, vcount))
#its common to center in the largest component
yeast.gc <- decompose.graph(yeast)[[1]]
average.path.length(yeast.gc)
diameter(yeast.gc)
transitivity(yeast.gc)

# k-vertex-connected if the removal of any subset of vertices of cardinality < k
#leaves a subgraph that is connected
#vertex connectivity is the largest integer such that
#G is k-vertex- connected.
vertex.connectivity(yeast.gc)
edge.connectivity(yeast.gc)
#vertex connectivity is bounded above by edge connectivity, 
#which is bounded above by min degree 

#If the removal of a particular set of vertices (edges) disconnects the graph,
# that set is called a vertex-cut (edge-cut).
# A single vertex that disconnects the graph is called a cut vertex, or articulation point.
# Identification of such vertices can provide a sense of where a network is vulnerable
yeast.cut.vertices <- articulation.points(yeast.gc)
length(yeast.cut.vertices)/vcount(yeast)
#[1] 0.1337409

#for directed nets
is.connected(aidsblog, mode=c("weak"))
is.connected(aidsblog, mode=c("strong"))
#maximal (weakly or strongly) connected components of a graph
aidsblog.scc <- clusters(aidsblog, mode=c("strong"))#same output than components()
table(aidsblog.scc$csize)

#assortativity=selective linking among vertices~correlation coefficients
#assortative characteristics can be categorical, ordinal, or continuous
assortativity.nominal(yeast, (V(yeast)$Class=="P")+1, directed=FALSE)
#One common use of the assortativity coefficient is summarizing degree–degree correlation of adjacent v
assortativity.degree(yeast)

######################PARTITIONING########################
library(ape)

#hierarchical clustering
kc <- fastgreedy.community(karate)
length(kc)
sizes(kc)
membership(kc)
par(mfrow=c(1,2))
#color clusters in graph
plot(kc, karate)
#dendrogram with clusters highlighted 
dendPlot(kc, mode="phylo")

#spectral clustering
k.lap <- graph.laplacian(karate)
#The graph Laplacian of a graph G, with adjacency matrix A, is a matrix L = D − A, 
#where D is a diagonal matrix with the degree sequence 
eig.anal <- eigen(k.lap)
plot(eig.anal$values, col="blue", ylab="Eigenvalues of Graph Laplacian")
#the smaller λ2 the more amenable the graph is to being separated into two subgraphs
f.vec <- eig.anal$vectors[, 33]#33=Fiedler vector:the 1st>0 (from the 34 decreasing λ)
f.colors <- ifelse(V(karate)$Faction == 1, "red", "cyan")
plot(f.vec, pch=16, xlab="Actor Number",
 ylab="Fiedler Vector Entry", col=f.colors)
abline(0, 0, lwd=2, col="lightgray")
#Fiedler suggested separating vertices according to their sign in the vector

#validation
func.class <- get.vertex.attribute(yeast.gc, "Class")
table(func.class)
yc <- fastgreedy.community(yeast.gc)
c.m <- membership(yc)
#confussion table
table(c.m, func.class, useNA=c("no"))