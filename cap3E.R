library(igraph)
library(sand)

data(Ecoli.data)
data(ppi.CC)
g=graph.adjacency(regDB.adj)

plot(g, layout=layout.circle)
#fea
plot(g, layout=layout.fruchterman.reingold)
#dos componentes aparentes
plot(g, layout=layout.kamada.kawai)
#varios componentes
plot(g, layout=layout.reingold.tilford)
#graph contains a cycle

#Decorating
plot(g, layout=layout.kamada.kawai,
	vertex.label=sapply(strsplit(V(g)$name,"_"),function(x) x[1]))

nbhds <- graph.neighborhood(g, order=1)
#neighborhood of each node
sapply(nbhds, vcount)
which(sapply(nbhds, vcount)>30)
plot(nbhds[[22]],
	vertex.label=sapply(strsplit(V(nbhds[[22]])$name,"_"),function(x) x[1]),
	layout=layout.kamada.kawai)

plot(ppi.CC, layout=layout.reingold.tilford(ppi.CC,circular=T))
#layout_nicely chooses automatically one