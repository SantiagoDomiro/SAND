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
 plane(state=s[1,],parms=p,odes=model,legend=T,xmax=2)
 abline(v=0.5,lty=2)
 abline(h=1.0368,lty=2)
dev.off()