library(sand)
library(ergm) # Will load package ’network’ as well.

data(lazega)
A <- as_adjacency_matrix(lazega)
v.attrs <- as_data_frame(lazega, what="vertices")
#transform data to proper format
lazega.s <- as.network(as.matrix(A),directed=FALSE)
set.vertex.attribute(lazega.s, "Office", v.attrs$Office)
set.vertex.attribute(lazega.s, "Practice",v.attrs$Practice)
set.vertex.attribute(lazega.s, "Gender",v.attrs$Gender)
set.vertex.attribute(lazega.s, "Seniority",v.attrs$Seniority)

###################################################fit models
set.seed(42)
#simpler model, depends on edge number
my.ergm.bern <- formula(lazega.s ~ edges)
simpler.fit <- ergm(my.ergm.bern)
#check deviance, not p-val coz chisq aint proven for this
anova(simpler.fit)
#Model 1: lazega.s ~ edges
#         Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
#NULL                       630     873.37                 
#Model 1:  1   274.58       629     598.78    < 2.2e-16 ***
#model considering kstars and triangles
my.ergm.kstar.tri <- formula(lazega.s ~ edges + kstar(2)
   + kstar(3) + triangle)
kstar.tri.fit <-ergm(my.ergm.kstar.tri)
#Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.
anova(kstar.tri.fit)    

#model with geometrically weighted degree count
my.ergm.gwd <- formula(lazega.s ~ edges
   + gwesp(1, fixed=TRUE))
gwd.fit <- ergm(my.ergm.gwd)
#Iteration 2 of at most 20:
#The log-likelihood improved by 7.237.
#..Iteration 13 of at most 20:
#The log-likelihood improved by 0.01073. <-stop at convergence
anova(gwd.fit)
#Model 1: lazega.s ~ edges + gwesp(1, fixed = TRUE)
#         Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
#NULL                       630     873.37                 
#Model 1:  2   352.39       628     520.97    < 2.2e-16 ***

#model with node characteristics, think on assortativity
lazega.ergm <- formula(lazega.s ~ edges
   + gwesp(log(3), fixed=TRUE)
   + nodemain("Seniority")#sum of attr(i) and attr(j) for all edges (i,j)
   + nodemain("Practice")#sum of attr(i) and attr(j) for all (i,j)
   + match("Practice") #homophily
   + match("Gender") #homophily
   + match("Office")) #homophily
#https://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html#term_nodemain_3
lazega.ergm.fit <- ergm(lazega.ergm)
anova(lazega.ergm.fit)
## Model 1: lazega.s ~ edges + gwesp(log(3), fixed = TRUE) + 
##     nodemain("Seniority") + nodemain("Practice") + 
##     match("Practice") + match("Gender") + 
##     match("Office")
##          Df Deviance Resid. Df Resid. Dev Pr(>|Chisq|)    
## NULL                       630     873.37                 
## Model 1:  7   413.74       623     459.63    < 2.2e-16 ***
#best fit

# the estimated coefficients in this analysis may be interpreted
# as a conditional log-odds ratio for cooperation between lawyers
summary(lazega.ergm.fit)
## Monte Carlo MLE Results:
##                              Estimate Std. Error MCMC %      
## edges                        -7.00655    0.67114      0   
## gwesp.fixed.1.09861228866811  0.59166    0.08554      0     
## nodecov.Seniority             0.02456    0.00620      0     
## nodecov.Practice              0.39455    0.10218      0    
## nodematch.Practice            0.76966    0.19060      0     
## nodematch.Gender              0.73767    0.24362      0     
## nodematch.Office              1.16439    0.18753      0     
##
##                               z value     Pr(>|z|)
## edges                         -10.440    < 1e-04 ***
## gwesp.fixed.1.09861228866811    6.917    < 1e-04 ***
## nodecov.Seniority               3.962    < 1e-04 ***
## nodecov.Practice                3.861    0.000113 ***
## nodematch.Practice              4.038    < 1e-04 ***
## nodematch.Gender                3.028    0.002463 **
## nodematch.Office                6.209    < 1e-04 ***
## 
## AIC: 473.6    BIC: 504.7    (Smaller is better.)

#check exp(Estimate)
coefs=summary(lazega.ergm.fit)$coefficients
cbind(rownames(coefs),exp(coefs$Estimate))
#"edges"                        "0.000906712982935729"
#"gwesp.fixed.1.09861228866811" "1.80290539705403"   
#"nodecov.Seniority"            "1.02492736852072"    
#"nodecov.Practice"             "1.48526392579601"   
#"nodematch.Practice"           "2.15948510062075"   2xodds of coperation 
#"nodematch.Gender"             "2.09445709065698"   2xodds of coperation 
#"nodematch.Office"             "3.18372373696905"   3xodds of coperation
#assortativity is driven by office

####################################################godness of fit
gof.lazega.ergm <- gof(lazega.ergm.fit)
plot(gof.lazega.ergm)