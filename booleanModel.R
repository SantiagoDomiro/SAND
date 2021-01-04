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
	add.expr = abline(v=which(basins[[x]][,20]==0),
	col="red",lwd=2),margins=c(1,8),
	labRow=gsub("initialState.","",rownames(temp[[x]]))))
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
#save colors for coming plots
colors=plotStateGraph(attrAND,layout=layout.svd,
	colorsAlpha=c(.3,.9,1,1),plotIt=F)
colors=unique(V(colors)$color)

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
get_basins<-function(model){
	sim<- markovSimulation(model)#get attractors
	g=plotPBNTransitions(sim,plotIt=F)#igraph
	basins=components(g)$membership#get basins
	#format 
	basins=sapply(1:8,function(x) names(which(basins==x)))
	basins=lapply(1:8,function(x) 
		sapply(strsplit(basins[,x],'*'),as.numeric))
	return(list(sim=sim,basins=basins,g=g))
}
zebProb_r=get_basins(zebProb)

#compare to non-probabilistic basins
sapply(zebProb_r$basins,ncol)
#[1] 64 64 64 64 64 64 64 64
sapply(1:8,function(x) sum(apply(zebProb_r$basins[[x]],2,
	paste,collapse="")%in%apply(temp[[x]],2,paste,collapse="")))
#[1] 64 64 64 64 64 64 64 64#same states per basin

#are attractors the same????
V(zebProb_r$g)$type=substr(V(zebProb_r$g)$name,5,6)%in%c("10","01","11")
#plot per attractor
pdf("Downloads/zeb1Proba_attr.pdf")
par(mar=c(1,1,1,5))
lapply(1:8,function(x) {p=induced_subgraph(zebProb_r$g,
	neighbors(zebProb_r$g,
	paste(temp[[x]][,basins[[x]]$transitionsToAttractor==0],
	collapse="")));
	plot(p,vertex.color=c("cornflowerblue","tomato")[
		as.factor(V(p)$type)],edge.width=as.numeric(E(p)$label)*10,
		vertex.label.color="black",edge.label.color="black",
		layout=layout.circle)})
dev.off()
#plot all together
col=colors[grep("4D",colors)]
V(zebProb_r$g)$membership=col[as.factor(components(zebProb_r$g)$
	membership)]
E(zebProb_r$g)$membership=sapply(get.edgelist(zebProb_r$g)[,1],
	function(x) V(zebProb_r$g)$membership[V(zebProb_r$g)$name==x])
png("Downloads/zebProb_basins.png",width=800,height=800)
par(mar=c(1,1,1,11))
plot(zebProb_r$g,vertex.size=3,vertex.label=NA,edge.label=NA,
 edge.width=as.numeric(E(zebProb_r$g)$label)*7,
 vertex.color=V(zebProb_r$g)$membership,
 edge.arrow.size=0.4,vertex.frame.color=NA,
 edge.color=col[as.factor(E(zebProb_r$g)$membership)],
 layout=layout.fruchterman.reingold)
dev.off()
#####################################model with 2 DICER copies
zeb2Dicer=loadNetwork("Downloads/ZEB12Dicer.txt")#2 copies of DICER
zeb2Dicer_r=get_basins(zeb2Dicer)
knockedDicer <- fixGenes(zeb2Dicer, "Dicer_1", 0)#1 mutated
knockedDicer_r=get_basins(knockedDicer)

#compare 1 hit vs non-hit model
sapply(zeb2Dicer_r$basins,ncol)
#[1] 128 128 128 128 128 128 128 128
sapply(knockedDicer_r$basins,ncol)
#[1] 64 64 64 64 64 64 64 64
sapply(1:8,function(x) sum(apply(zeb2Dicer_r$basins[[x]],2,paste,
	collapse="")%in%apply(knockedDicer_r$basins[[x]],2,paste,
	collapse="")))

#compare 1 hit vs 1 copy model
sapply(1:8,function(x) sum(apply(knockedDicer_r$basins[[x]][2:10,],2,
	paste,collapse="")%in%apply(temp[[x]],2,paste,collapse="")))
#[1] 64 64 64 64 64 64 64 64
#plot basins
#knockedDicer_r$basins=lapply(1:8,function(x) 
#	knockedDicer_r$basins[[x]][,
#	order(match(apply(knockedDicer_r$basins[[x]][2:10,],2,paste,
#		collapse=""),apply(temp[[x]],2,paste,collapse="")))])
#lapply(1:8,function(x) heatmap(knockedDicer_r$basins[[x]],Rowv=NA,
#	scale="n",col=c("white","black"),Colv=NA,labCol=NA,margins=c(1,8),
#	labRow=knockedDicer_r$sim$genes))

#plot attractors
V(knockedDicer_r$g)$type=substr(V(knockedDicer_r$g)$name,6,
	7)%in%c("10","01","11")#update to match EMT & resist positions
pdf("Downloads/knockedDicer_attr.pdf")
par(mar=c(1,1,1,5))
lapply(1:8,function(x) {p=induced_subgraph(knockedDicer_r$g,
	neighbors(knockedDicer_r$g,
	paste(0,paste(temp[[x]][,basins[[x]]$transitionsToAttractor==0],
	collapse=""),sep="")));
	plot(p,vertex.color=c("cornflowerblue","tomato")[
		as.factor(V(p)$type)],edge.width=as.numeric(E(p)$label)*10,
		vertex.label.color="black",edge.label.color="black",
		layout=layout.circle)})
dev.off()


#possible follwing steps
#1 strong vs weak basins of atraction (for assync)
#2 entropy of the net
#propagation of info
#net into discrete time bi-linear system 
# stability of attractors
#perturbations

