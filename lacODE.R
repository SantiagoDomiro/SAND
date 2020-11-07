#modified from http://tbb.bio.uu.nl/rdb/grindR/operon.R
source("Documents/scripts/grind.R")

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
continue(state=s[1,],parms=p,odes=model,x="L",y="M",xmax=2,ymax=1.3,color="black")
dev.off()
