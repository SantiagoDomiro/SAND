library(sand)

set.seed(42)
g.er <- erdos.renyi.game(100, 0.02)#simulate classical random graph
#100 nodes and prob of edge=0.02
plot(g.er, layout=layout.circle, vertex.label=NA)

#####################################properties of random graphs
#graph is not necesarily connected
is.connected(g.er)
# [1] FALSE

#but it has a giant component
table(sapply(decompose.graph(g.er), vcount))
# 1  2  3  4 71 
#15  2  2  1  1 

#mean degree is close to expected (100 − 1) × 0.02 = 1.98
mean(degree(g.er))
#[1] 1.9
#distribution is homogeneous
hist(degree(g.er), col="lightblue",xlab="Degree",
 ylab="Frequency", main="")

#short shortest paths
average.path.length(g.er)
# [1] 5.276511
diameter(g.er)
# [1] 14

#low clustering.
transitivity(g.er)
#[1] 0.01639344
#diameter scales like log V & clustering behaves like V**−1
########################################null models conserving degree
#sample random graphs with a fixed degree sequence
degs <- c(2,2,2,2,3,3,3,3)
g1 <- degree.sequence.game(degs, method="vl")
g2 <- degree.sequence.game(degs, method="vl")
par(mfrow=c(1,2))
plot(g1, vertex.label=NA)
plot(g2, vertex.label=NA)

#gs are not relabellings of each other but share ecount
graph.isomorphic(g1, g2)
# [1] FALSE
c(ecount(g1), ecount(g2))
#[1] 10 10


#sample random graphs with a real degree sequence
data(yeast)
degs <- degree(yeast)
fake.yeast <- degree.sequence.game(degs,method=c("vl"))
#all(): Given a set of logical vectors, are all of the values true?
all(degree(yeast) == degree(fake.yeast))
#[1] TRUE


#non-degree features vary
diameter(yeast)
#[1] 15
diameter(fake.yeast)#half as large
#[1] 8
transitivity(yeast)
#[1] 0.4686178
transitivity(fake.yeast)#quite lower
#[1] 0.04026804


######################################## Small-World Models
g.ws <- sample_smallworld(1, 25, 5, 0.05)#(dim,vertex,neighbors,p)
plot(g.ws, layout=layout_in_circle, vertex.label=NA)

#lattice vs rewiring
g.lat100 <- sample_smallworld(1, 100, 5, 0)
transitivity(g.lat100)
#[1] 0.6666667
diameter(g.lat100)
#[1] 10
mean_distance(g.lat100)
#[1] 5.454545

g.ws100 <- sample_smallworld(1, 100, 5, 0.05)
transitivity(g.ws100)
#[1] 0.4974675 a little bit lower
diameter(g.ws100)
#[1] 5
mean_distance(g.ws100)# average shortest paths between all pairs of v 
#[1] 2.725051 a little bit higher

steps <- seq(-4, -0.5, 0.1)
len <- length(steps)
ntrials <- 100
#mean cl and dist from 100 small worlds per p=100**step[i]
alter<-sapply(1:len,function(i) rowMeans(sapply(1:ntrials,function(j)
 {g=sample_smallworld(1, 1000, 10, 10**steps[i]);
 	return(cbind(transitivity(g),mean_distance(g)))})))

plot(steps, alter[1,]/max(alter[1,]), ylim=c(0, 1), lwd=3, type="l",
 col="blue", xlab=expression(log[10](p)),
 ylab="Clustering and Average Path Length")
lines(steps, alter[2,]/max(alter[2,]), lwd=3, col="red")
#cl decreases slowly first and fast later
#mean distance does the opposite

########################### Preferential Attachment Models
set.seed(42)
g.ba <- sample_pa(100, directed=FALSE)#100 added nodes
#the larger t, the larger difference among regular nodes & hubs
rowMeans(sapply(1:10,function(x) 
	summary(hub.score(sample_pa(10,directed=F))$vector)))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.1743031 0.2606717 0.3384106 0.4480632 0.5882570 1.0000000 
rowMeans(sapply(1:10,function(x) 
	summary(hub.score(sample_pa(100,directed=F))$vector)))
#       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
#0.000652187 0.015384800 0.041521715 0.092000482 0.129205265 1.000000000 
rowMeans(sapply(1:10,function(x) 
	summary(hub.score(sample_pa(1000,directed=F))$vector)))
#        Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
#5.190090e-08 1.758885e-04 1.297652e-03 1.352760e-02 7.855621e-03 1.000000e+00 

#properties 
par(mfrow=c(1,2))
plot(g.ba, layout=layout_in_circle, vertex.label=NA)
hist(degree(g.ba), col="lightblue",
 xlab="Degree", ylab="Frequency", main="")
summary(degree(g.ba))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    1.00    1.00    1.98    2.00    9.00 
#75 vertex have degree<2
mean_distance(g.ba)#short distances similar to lattices
#[1] 5.815556
diameter(g.ba)
#[1] 12
transitivity(g.ba)#low transitivity
#[1] 0
#Scale-free networks with degree exponents 2<γ<3 are ultra-small
#with the average path length following l ~ loglogN
#Barabási–Albert model does not have an inherent modularity,
# so C(k) is independent of k 

########################### Comunnities significance
# data(karate)
nv <- vcount(karate)
ne <- ecount(karate)
degs <- degree(karate)
ntrials <- 1000
#random graphs with the same vcount & ecount
cl_distri_rnd <- sapply(1:ntrials,function(x) {
 g.rg <- sample_gnm(nv, ne)
 c.rg <- cluster_fast_greedy(g.rg)
return(length(c.rg))})
#random with the same degree distri
cl_distri_d <- sapply(1:ntrials,function(x) {
 g.grg <- sample_degseq(degs, method="vl")
 c.grg <- cluster_fast_greedy(g.grg)
return(length(c.grg))})
#get distri of communities in random graphs
rslts <- c(cl_distri_rnd,cl_distri_d)
indx <- c(rep(0, ntrials), rep(1, ntrials))
counts <- table(indx, rslts)/ntrials
barplot(counts, beside=TRUE, col=c("blue", "red"),
 xlab="Number of Communities",
 ylab="Relative Frequency",
 legend=c("Fixed Size", "Fixed Degree Sequence"))

############################am I smallW????????
library(igraphdata)
data(macaque)
#implement cc for directed nets p.82
clust_coef_dir <- function(graph) {
   A <- as.matrix(as_adjacency_matrix(graph))
   S <- A + t(A)
   deg <- degree(graph, mode=c("total"))
   num <- diag(S %*% S %*% S)
   denom <- diag(A %*% A)
   denom <- 2 * (deg * (deg - 1) - 2 * denom)
   cl <- mean(num/denom)
return(cl)}

#cc from random similar graphs
nv <- vcount(macaque)
ne <- ecount(macaque)
rnd_propes <- t(sapply(1:ntrials,function(x) {
 g.rg <- sample_gnm(nv, ne, directed=TRUE)
return(cbind(clust_coef_dir(g.rg),mean_distance(g.rg)))}))
summary(rnd_propes)
#       V1               V2       
# Min.   :0.2169   Min.   :1.808  
# 1st Qu.:0.2303   1st Qu.:1.827  
# Median :0.2339   Median :1.833  
# Mean   :0.2342   Mean   :1.833  
# 3rd Qu.:0.2379   3rd Qu.:1.838  
# Max.   :0.2525   Max.   :1.862  

#not so similar
clust_coef_dir(macaque)
#[1] 0.5501073 larger than expected
mean_distance(macaque)
#[1] 2.148485 #also larger
#evidence for small-world behavior in this network is not clear
#better null model??????????