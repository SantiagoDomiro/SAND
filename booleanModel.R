library(BoolNet)
library(igraph)
#get model PMC4614883 & Rupaimoole
zebAND<-loadNetwork("Downloads/ZEB1AND.txt")#AND regulators
zebOR<-loadNetwork("Downloads/ZEB1.txt")#OR regulators
#as in https://cran.r-project.org/web/packages/BoolNet/vignettes/BoolNet_package_vignette.Snw.pdf
##############################################33333#sync analysis
#getPathToAttractor(zeb,rep(1,9))#states to attractor

#get all attractors
attrAND <- getAttractors(zebAND)
#print(attrAND, activeOnly=TRUE)
#recover all states in cancer basins
basins=lapply(1:8,function(x) 
	as.data.frame(getBasinOfAttraction(attrAND,x)))
#plot them
temp=lapply(basins,function(x) t(as.matrix(x[,1:9])))#columns with states
pdf("zeb1Basins.pdf")
lapply(1:8,function(x) heatmap(temp[[x]],Rowv=NA,scale="n",
	col=c("white","black"),Colv=NA,labCol=NA,
	add.expr = abline(v=which(cancer[[x]][,20]==0),
	col="red",lwd=2)))
dev.off()
#check hypoxia removal
getPathToAttractor(zebAND,rep(0,9))
#  Dicer miR200 miR205 ZEB1 EMT radioresistance hypoxia IR MeS
#     0      0      0    0   0               0       0  0   0
#     1      0      0    0   0               0       0  0   0
#     1      1      1    0   0               0       0  0   0

#get all attractors
attrOR <- getAttractors(zebOR)
print(attrOR, activeOnly=TRUE)

#getTransitionTable(attr) state transitions
#getBasinOfAttraction(attr, 2) #states in a basin
png("Downloads/ZEB1attr.png")
plotStateGraph(attrAND,layout=layout.svd,colorsAlpha=c(.3,.9,1,1))
dev.off()


#############################################async
attrAND.1 <- getAttractors(zebAND,
 type="asynchronous",
 method="random",
 startStates=500,
 avoidSelfLoops=FALSE)
print(attrAND.1, activeOnly=TRUE)
#same as in synchronous
#####################################probabilistic model
zebProb=loadNetwork("Downloads/ZEB1.txt")#Each functions has a probability, probs per var sum up to 1
#get attractors
sim <- markovSimulation(zebProb)
g=plotPBNTransitions(sim,plotIt=F)#igraph
#get basins
basinsProb=components(g)$membership
table(basinsProb)
# 1  2  3  4  5  6  7  8 
#64 64 64 64 64 64 64 64 #same number of attractors & size of basins

#format to plot
basinsProb=sapply(1:8,function(x) names(which(basinsProb==x)))
basinsProb=lapply(1:8,function(x) 
	sapply(strsplit(basinsProb[,x],'*'),as.numeric))
#compare to non-probabilistic basins
sapply(1:8,function(x) sum(apply(basinsProb[[x]],2,paste,collapse="")
	%in%apply(temp[[x]],2,paste,collapse="")))
#[1] 64 64 64 64 64 64 64 64#same states per basin
#pdf("Downloads/zeb1BasinsProba.pdf")
#lapply(1:8,function(x) heatmap(basins[[x]],Rowv=NA,scale="n",
#	col=c("white","black"),Colv=NA,labCol=NA))
#dev.off()

#are attractors the same????
V(g)$type=substr(V(g)$name,5,6)%in%c("10","01","11")
#plot per attractor
pdf("Downloads/zeb1BasinsProba.pdf")
lapply(1:8,function(x) {p=induced_subgraph(g,neighbors(g,
	paste(basins[[x]][basins[[x]]$transitionsToAttractor==0,1:9],
	collapse="")));
	plot(p,vertex.color=c("cornflowerblue","tomato")[
		as.factor(V(p)$type)],edge.width=as.numeric(E(p)$label)*10,
		vertex.label.color="black",edge.label.color="black")})
dev.off()

#possible follwing steps
#1 strong vs weak basins of atraction (for assync)
#2 entropy of the net
#propagation of info
#net into discrete time bi-linear system 
# stability of attractors
#perturbations

