source("Documents/scripts/grind.R")

#modified from http://tbb.bio.uu.nl/rdb/grindR/operon.R
model <- function(t, state, parms) {
   with(as.list(c(state,parms)), { 
   R = 1/(1+A^n)               # Repressor
   dM = c0 + c*(1-R) - d*M     # mRNA
   dA = M*L - delta*A - (u*M*A)/(h+A)  # Allolactose
   return(list(c(dA, dM)))  
})}

p <- c(L=1,c=1,d=1,u=1,c0=0.05,h=2,n=5,delta=0.2)#params
s <- cbind(A=c(seq(0,1,length=4),4,0),M=c(seq(0,1,length=4),0,4))#initial states

#time plot
png("Downloads/time.png",width=800)
 par(mfcol=c(2,3))
 apply(s,1,function(x) round(run(tmax=40,
 								tstep=1,
							 	odes=model,
							 	state=x,
							 	parms=p),4))
#    [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
#A 0.2271 0.2274 2.3717 2.3717 2.3717 2.3717
#M 0.0506 0.0506 1.0368 1.0368 1.0369 1.0368
dev.off()

#add noise on L concentration
png("Downloads/time+noise.png")
par(mfrow=c(2,2))
sapply(seq(0.5,2,0.5),function(l) {
	p["L"]=l;
	run(tmax=40,tstep=1,odes=model,state=s[4,],parms=p,
		after="parms[\"c\"]<-rnorm(1,mean=1,sd=0.2)",
		main=paste("L =",p["L"]))})
dev.off()
#        [,1]     [,2]     [,3]     [,4]
#A 0.11213256 2.444345 4.197148 6.552934
#M 0.05001575 1.210919 0.969043 1.031852

#phase plane
png("Downloads/phase.png")
 plane(state=s[1,],parms=p,odes=model,legend=T,xmax=3,ymax=1.2,
  vector=T,grid=14)# find steady states
 low <- newton(s[1,],plot=T)#s[1,]=0.0000000 0.0000000
 #         A          M 
 #0.22721312 0.05060519 
 #Stable point, eigenvalues:  -1.015041 -0.2053624 
 mid <- newton(s[4,],plot=T)#s[4,]=1.0000000 1.0000000
 #        A         M 
 #0.6907057 0.1858486 
 #Unstable point, eigenvalues:  -1.504184 0.2528436  
 high <- newton(s[5,],plot=T)#s[5,]=4.0000000 0.0000000
 #      A       M 
 #2.37172 1.03685 
 #Stable point, eigenvalues:  -1.01765 -0.2908533 
dev.off()

#################varying suggested parameter
p["L"]<-0.5
png("Downloads/phaseAlt.png")
plane(state=s[1,],parms=p,odes=model,legend=T,xmax=3,ymax=1.2,
 vector=T,grid=14,main="L = 0.5")# find steady states
low <- newton(s[1,],plot=T)#s[1,]=0.0000000 0.0000000
#         A          M 
#0.11180350 0.05001747 
#Stable point, eigenvalues:  -1.000449 -0.2219819 
mid <- newton(s[4,],plot=T)#s[4,]=1.0000000 1.0000000
#         A          M 
#0.11180350 0.05001747 
#Stable point, eigenvalues:  -1.000449 -0.2219819 
high <- newton(s[5,],plot=T)#s[5,]=4.0000000 0.0000000
#         A          M 
#-0.6893587 -0.1343818 #those value make no sense
#Unstable point, eigenvalues:  -1.8833 0.8397599 
dev.off()

png("Downloads/bifurcationM.png")
continue(state=s[1,],parms=p,odes=model,x="L",y="M",
	xmax=2,ymax=1.3,color="black")
dev.off()

######################################check repressor function
repre<-function(A,n){1/(1+A^n)}
repreData=data.frame(do.call(rbind,lapply(seq(0.1,3,length=20),
	function(A) t(sapply(1:10,function(n) 
		cbind(A,n,repre(A,n)))))))
colnames(repreData)=c("A","n","R")
png("Downloads/hillR.png")
 ggplot(repreData,aes(x=A,y=R,col=as.character(n)))+geom_point()+
 geom_line()+labs(color="n")
dev.off()

modelExtended <- function(t, state, parms) {
   with(as.list(c(state,parms)), { 
   dR = p-r*R-1/(1+A^n)              # Repressor
   dM = c0 + c*(1-R) - d*M     # mRNA
   dA = M*L - delta*A - (u*M*A)/(h+A)  # Allolactose
   return(list(c(dR, dM, dA)))  
})}
p <- c(p,c(p=1,r=0.4))#params
s <- c(R=1,M=0,A=0)#initial states

modelAlt <- function(t, state, parms) {
   with(as.list(c(state,parms)), { 
   dR = R-1/(1+A^n)              # Repressor
   dM = c0 + c*(1-R) - d*M     # mRNA
   dA = M*L - delta*A - (u*M*A)/(h+A)  # Allolactose
   return(list(c(dR, dM, dA)))  
})}
png("Downloads/Extended.png")
par(mfrow=c(2,3))
run(odes=modelExtended,state=s,parms=p,tstep=1,tmax=150)
#        R         M         A 
#0.8041728 0.2462669 0.8616264 
#plane(state=s,parms=p,odes=modelExtended,vector=T,grid=10,x=3,y=2,xmax=3)# find steady states
plane(state=s,parms=p,odes=modelExtended,vector=T,grid=10,x=3,y=1,ymax=3,xmax=3)
newton(s,plot=T,odes=modelExtended,x=3,y=1)
#        R         M         A 
#0.8035953 0.2464047 0.8611952 
#Stable point, eigenvalues:  -1.573951+0i -0.0431238+0.7916238i -0.0431238-0.7916238i 
continue(state=s,parms=p,odes=modelExtended,x="L",y="M")
#Starting at L = 1 with:
#        R         M         A 
#0.8035953 0.2464047 0.8611952 
#Bifurcation at L = 1.2 
run(odes=modelAlt,state=s,parms=p,tstep=1,tmax=40)
#        R         M         A 
#-4.001359 -1.093514 -1.000000 
#plane(state=s,parms=p,odes=modelExtended,vector=T,grid=10,x=3,y=2,xmax=3)# find steady states
plane(state=s,parms=p,odes=modelAlt,vector=T,grid=10,x=3,y=1,xmax=3)# find steady states
newton(s,plot=T,odes=modelAlt,x=3,y=1)
#         R          M          A 
#0.99939481 0.05060519 0.22721312 
#Unstable point, eigenvalues:  -1.007563 0.9950712 -0.2079112 
newton(c(R=0.2,M=3,A=3),plot=T,odes=modelAlt,x=3,y=1)
#         R          M          A 
#0.01315024 1.03684976 2.37172027 
#Unstable point, eigenvalues:  -1.008895 0.9951882 -0.2947955
continue(state=s,parms=p,odes=modelAlt,x="L",y="M")
#Starting at L = 1 with:
#         R          M          A 
#0.99939481 0.05060519 0.22721312 
#Turning point point at L = 1.510469 
#Turning point point at L = 0.6846875 
dev.off()

