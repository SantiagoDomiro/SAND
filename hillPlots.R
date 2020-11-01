library(ggplot2)
library(gridExtra)

hill_activator<-function(beta_max,activator,K,n){
return((beta_max*activator**n)/(K**n+activator**n))}

varia_n=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x) 
	do.call(rbind,lapply(c(seq(1,10,2),25,50,100,500),function(y) 
		cbind(hill_activator(1,x,0.5,y),x,y))))))
colnames(varia_n)=c("fx","x","n")
varia_n$n=factor(varia_n$n)

plotn=ggplot(varia_n,aes(x=x,y=fx,col=n))+geom_point()+
theme(text=element_text(size=18))+ggtitle("beta_max=1,K=0.5")

varia_k=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x)
 do.call(rbind,lapply(seq(0.25,0.75,length=10),function(k) 
 	cbind(hill_activator(1,x,k,5),x,k))))))
colnames(varia_k)=c("fx","x","k")

plotk=ggplot(varia_k,aes(x=x,y=fx,col=k))+geom_point()+
theme(text=element_text(size=18))+ggtitle("beta_max=1,n=5")

varia_beta=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x) 
	do.call(rbind,lapply(seq(0.125,1,0.125),function(beta) 
		cbind(hill_activator(beta,x,0.5,5),x,beta))))))
colnames(varia_beta)=c("fx","x","beta")

plotb=ggplot(varia_beta,aes(x=x,y=fx,col=beta))+geom_point()+
theme(text=element_text(size=18))+ggtitle("k=0.5,n=5")

png("Desktop/activador.png",width=1000)
grid.arrange(plotn,plotk,plotb,ncol=3)
dev.off()

hill_represor<-function(beta_max,represor,K,n){
	return(beta_max/(1+(represor/K)**n))}

varia_n=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x) 
	do.call(rbind,lapply(c(seq(1,10,2),25,50,100,500),function(y) 
		cbind(hill_represor(1,x,0.5,y),x,y))))))
colnames(varia_n)=c("fx","x","n")
varia_n$n=factor(varia_n$n)
plotn=ggplot(varia_n,aes(x=x,y=fx,col=n))+geom_point()+
theme(text=element_text(size=18))+ggtitle("beta_max=1,K=0.5")

varia_k=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x)
 do.call(rbind,lapply(seq(0.25,0.75,length=10),function(k) 
 	cbind(hill_represor(1,x,k,5),x,k))))))
colnames(varia_k)=c("fx","x","k")
plotk=ggplot(varia_k,aes(x=x,y=fx,col=k))+geom_point()+
theme(text=element_text(size=18))+ggtitle("beta_max=1,n=5")

varia_beta=data.frame(do.call(rbind,lapply(seq(0,1.5,length=100),function(x) 
	do.call(rbind,lapply(seq(0.125,1,0.125),function(beta) 
		cbind(hill_represor(beta,x,0.5,5),x,beta))))))
colnames(varia_beta)=c("fx","x","beta")
plotb=ggplot(varia_beta,aes(x=x,y=fx,col=beta))+geom_point()+
theme(text=element_text(size=18))+ggtitle("k=0.5,n=5")

png("Desktop/represor.png",width=1000)
grid.arrange(plotn,plotk,plotb,ncol=3)
dev.off()
######################
protein=function(beta,alpha,time){return(beta/alpha*(1-expm1(1)**(-alpha*time)))}
prote_alpha=data.frame(do.call(rbind,lapply(1:100,function(t) 
	do.call(rbind,lapply(seq(0.1,1,length=10),function(alpha)
	 cbind(t,protein(1,alpha,t),alpha))))))
prote_alpha$t=as.numeric(as.character(prote_alpha$t))
colnames(prote_alpha)[2]="ft"
prote_alpha$ft=as.numeric(as.character(prote_alpha$ft))
pa=ggplot(prote_alpha,aes(x=t,y=ft,col=alpha))+geom_point()+
ggtitle("beta=1")+theme(text=element_text(size=18))

prote_beta=data.frame(do.call(rbind,lapply(1:100,function(t) 
	do.call(rbind,lapply(seq(0.1,1,length=10),function(beta)
	 cbind(t,protein(beta,0.5,t),beta))))))
prote_beta$t=as.numeric(as.character(prote_beta$t))
colnames(prote_beta)[2]="ft"
prote_beta$ft=as.numeric(as.character(prote_beta$ft))
pb=ggplot(prote_beta,aes(x=t,y=ft,col=beta))+geom_point()+
ggtitle("alpha=0.5")+theme(text=element_text(size=18))
png("Desktop/protein.png",width=800)
grid.arrange(pa,pb,ncol=2)
dev.off()

png("Desktop/t_small.png")
ggplot(prote_alpha,aes(x=t,y=ft,group=alpha,col=alpha))+geom_point()+xlim(0,10)+
ylim(0,5)+geom_smooth(method='lm',formula=y~x)
dev.off()