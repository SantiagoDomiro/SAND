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
print(attrAND, activeOnly=TRUE)
#recover all states in cancer basins
cancer=lapply(c(4,6,8),function(x) 
	as.data.frame(getBasinOfAttraction(attrAND,x)))
temp=lapply(cancer,function(x) t(as.matrix(x[,1:9])))#columns with states
heatmap(temp[[2]],Rowv=NA,scale="n",col=c("white","black"),
	Colv=NA,labCol=NA,add.expr = abline(v=which(cancer[[2]][,20]==0),
		col="red",lwd=2))

#get all attractors
attrOR <- getAttractors(zebOR)
print(attrOR, activeOnly=TRUE)

#getTransitionTable(attr) state transitions
#getBasinOfAttraction(attr, 2) #states in a basin
png("Downloads/ZEB1basin.png")
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
#####################################check DE
i=rbind(i,myannot[myannot$hgnc_symbol=="CDH1",])
DE.i=lapply(DE.genes,function(x) x[rownames(x)%in%i$ensembl_gene_id,])
DE.i=lapply(DE.i,function(x) x[order(match(rownames(x),i$ensembl_gene_id)),])
DE.i=do.call(cbind,lapply(DE.i,function(x) x$logFC))
rownames(DE.i)=i$hgnc_symbol
colnames(DE.i)=gsub("_Normal","",colnames(DE.i))

#possible follwing steps
#1 strong vs weak basins of atraction (for assync)
#2 entropy of the net
#propagation of info
#net into discrete time bi-linear system 
# stability of attractors
#perturbations

