#Code Samples N given W matrix

# #####################
# ####Simulate Data####
# #####################
# 
# setwd("C:/Users/Amanda/Dropbox/A_Ellis/Code")
# source("Sim_Data.R")                           #Source Data Simulation function
# 
# t=3
# p=.5
# N=30
# lambda=2
# alpha.match=1
# beta.match=1
# alpha.non.match=3
# beta.non.match=3
# 
# W<-sim.data.M0(t,p,N,alpha.n,beta.n,alpha.match,beta.match,alpha.non.match,beta.non.match,'W')  #Simulates Data

###########################################
#### Manipulate Data to Use in Sampler ####
###########################################

sampler.M0(W,interations=100000,burn.in=10000,CI=.95){
n.occasions<-length(W[1,])            #Number of capture Occasions
n.obs.ind<-length(W[,1])              #Number of observed individuals

captures<-rep(NA,length=n.occasions)  #Calculate the number of captures per capture occasion
for(i in 1:n.occasions){
  captures[i]<-sum(W[,i])
}

total.captures<-sum(captures)

first<-rep(NA,length=n.obs.ind) 
W.no.first<-W
for(i in 1:n.obs.ind){
  first[i]<-min(which(W[i,]==1))        #Calculate when individual first observed
  W.no.first[i,first[i]]<-0             #Sets the first occasion to 0
}

recaptures<-rep(NA,length=n.occasions)  #Calculate the number of recaptures per capture occasion
for(i in 1:n.occasions){
  recaptures[i]<-sum(W.no.first[,i])
}

################
#### MCMC  #####
################

#Define Chains
N.gibbs<-rep(NA,length=iterations+burn.in)
p.gibbs<-rep(NA,length=iterations+burn.in)

#intial values
N.gibbs[1]<-n.obs.ind
p.gibbs[1] <-.5

#Priors

alpha.p<-.5   #p has a beta prior
beta.p <-.5

alpha.N<-0    #N has negative binomial prior
beta.N <-0


for(i in 2:(iterations+burn.in)){                             #Gibbs sampler
  #Sample N  
  a=n.obs.ind+alpha.N
  b=(1+beta.N-((1-p.gibbs[i-1])^n.occasions))/((1-p.gibbs[i-1])^n.occasions)
  lambda<-rgamma(1,a,b)
  U<-rpois(1,lambda)
  N.gibbs[i]<- n.obs.ind + U 
  #N.gibbs[i]<-N
  
  #Sample p
  alpha.p.gibbs<-alpha.p+total.captures
  beta.p.gibbs<-n.occasions*N.gibbs[i]-total.captures+beta.p
  p.gibbs[i]<-rbeta(1,alpha.p.gibbs,beta.p.gibbs)
  #p.gibbs[i]<-p

}


N.gibbs <- N.gibbs[-(1:burn.in)]
#hist(N.gibbs, freq=FALSE, xlab="N", col="gray", border="white", main="")
posterior.mean.N<-mean(as.vector(N.gibbs),na.rm=TRUE)

p.gibbs <- p.gibbs[-(1:burn.in)]
#hist(p.gibbs, freq=FALSE, xlab="p", col="gray", border="white", main="")
posterior.mean.p<-mean(as.vector(p.gibbs),na.rm=TRUE)

prob=1-CI
lb=prob/2
ub=1-prob/2
CI.N<-quantile(N.gibbs,c(lb,ub),names=FALSE)
CI.p<-quantile(p.gibbs,c(lb,ub),names=FALSE)
return(list(posterior.man.N,posterior.mean.p,CI.N,CI.p))
}

