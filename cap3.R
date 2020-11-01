library(sand)

#5 × 5 × 5 lattice
g.l <- graph.lattice(c(5, 5, 5))
#blogs
data(aidsblog)
summary(aidsblog)

###################LAYOUTS##########################
#circular
igraph.options(vertex.size=3, vertex.label=NA,2,edge.arrow.size=0.5)
par(mfrow=c(1, 2))
plot(g.l, layout=layout.circle)
 title("5x5x5 Lattice")
plot(aidsblog, layout=layout.circle)
 title("Blog Network")

#spring-embedder methods
plot(g.l, layout=layout.fruchterman.reingold)
 title("5x5x5 Lattice")
plot(aidsblog, layout=layout.fruchterman.reingold)
 title("Blog Network")
#blog net looks different, seed?

#energy-placement methods
plot(g.l, layout=layout.kamada.kawai)
 title("5x5x5 Lattice")
plot(aidsblog, layout=layout.kamada.kawai)
 title("Blog Network")

#trees
g.tree <- graph.formula(1-+2,1-+3,1-+4,2-+5,2-+6,2-+7,3-+8,3-+9,4-+10)
igraph.options(vertex.size=30, edge.arrow.size=0.5,vertex.label=NULL)
par(mfrow=c(1, 3))
plot(g.tree, layout=layout.circle)
#radial
plot(g.tree, layout=layout.reingold.tilford(g.tree,circular=T))
#circular: whether to plot the tree in a circular fashion
plot(g.tree, layout=layout.reingold.tilford)

#bipartites
plot(g.bip, layout=-layout.bipartite(g.bip)[,2:1],
 vertex.size=30, vertex.shape=ifelse(V(g.bip)$type,
 "rectangle", "circle"),
 vertex.color=ifelse(V(g.bip)$type, "red", "cyan"))
#layout.bipartite gives u coordinates?

###################DECORATING GRAPHS###################
library(igraphdata)
data(karate)
#Reproducible layout
set.seed(42)
l <- layout.kamada.kawai(karate)
# Plot undecorated first.
igraph.options(vertex.size=10)
par(mfrow=c(1,1))
plot(karate, layout=l, vertex.label=V(karate))

# Now decorate, starting with labels.
V(karate)$label <- sub("Actor ", "", V(karate)$name)
# Two leaders get shapes different from club members.
V(karate)$shape <- "circle"
V(karate)[c("Mr Hi", "John A")]$shape <- "rectangle"
# Differentiate two factions by color.
V(karate)[Faction == 1]$color <- "red"
V(karate)[Faction == 2]$color <- "dodgerblue"
# Vertex area proportional to vertex strength
V(karate)$size <- 4*sqrt(graph.strength(karate))
#graph.strength:sum up edge weights of the adjacent edges for each vertex
V(karate)$size2 <- V(karate)$size * .5
# Weight edges by number of common activities
E(karate)$width <- E(karate)$weight
# Color edges by within/between faction.
F1 <- V(karate)[Faction==1]
F2 <- V(karate)[Faction==2]
E(karate)[ F1 %--% F1 ]$color <- "pink"
E(karate)[ F2 %--% F2 ]$color <- "lightblue"
E(karate)[ F1 %--% F2 ]$color <- "yellow"
#%--% is special operator used to select all edges between two sets of vertices
# Offset vertex labels for smaller points (default=0).
V(karate)$label.dist <- ifelse(V(karate)$size >= 10, 0, 0.75)
# Plot decorated graph, using same layout.
plot(karate, layout=l)

#egocentric viz
k.nbhds <- graph.neighborhood(karate, order=1)
#neighborhood of each node
sapply(k.nbhds, vcount)
#largest neighborhoods
k.1 <- k.nbhds[[1]]
k.34 <- k.nbhds[[34]]
par(mfrow=c(1,2))
 plot(k.1, vertex.label=NA,vertex.color=c("red", rep("lightblue", 16)))
 plot(k.34, vertex.label=NA,vertex.color=c(rep("lightblue", 17), "red"))

#supply parameters instead of setting them as attributes
library(sand)
data(lazega)
# Office location indicated by color.
colbar <- c("red", "dodgerblue", "goldenrod")
v.colors <- colbar[V(lazega)$Office]
# Type of practice indicated by vertex shape.
v.shapes <- c("circle", "square")[V(lazega)$Practice]
# Vertex size proportional to years with firm.
v.size <- 3.5*sqrt(V(lazega)$Years)
# Label vertices according to seniority.
v.label <- V(lazega)$Seniority
set.seed(42)
l <- layout.fruchterman.reingold(lazega)
plot(lazega, layout=l, vertex.color=v.colors,
 vertex.shape=v.shapes, vertex.size=v.size,
 vertex.label=v.label)

###################LARGE NETWORKS###################
summary(fblog)
party.names <- sort(unique(V(fblog)$PolParty))
party.names
set.seed(42)
l = layout.kamada.kawai(fblog)
party.nums.f <- as.factor(V(fblog)$PolParty)
party.nums <- as.numeric(party.nums.f)
plot(fblog, layout=l, vertex.label=NA,vertex.color=party.nums, vertex.size=3)

#meta-vertices
fblog.c <- contract.vertices(fblog, party.nums)#awesome
E(fblog.c)$weight <- 1
fblog.c <- simplify(fblog.c)#dump multiple edges
party.size <- as.vector(table(V(fblog)$PolParty))
plot(fblog.c, vertex.size=5*sqrt(party.size),vertex.label=party.names,
 vertex.color=V(fblog.c),edge.width=sqrt(E(fblog.c)$weight),
 vertex.label.dist=1.5, edge.arrow.size=0)
I