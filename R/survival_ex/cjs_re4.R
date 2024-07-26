rm(list=ls())
#data.dir<-"C:/users/mconroy/Dropbox/teaching/WILD8390/spring2015/exercises/wk13"


#setwd(data.dir)
# we will use JAGS instead of BUGS
require(jagsUI)
#simulate data under a CJS study with specified parameters
#see Kery and Schaub 2011 p 178

# VECTOR OF PHI AND P TO BE SPECIFIED LATER!!

#simulate CJS data for 1 group
simul.cjs<-function(phi,p,marked)
{
n.occasions<-length(p)+1
Phi<-matrix(phi,n.occasions-1,nrow=sum(marked),byrow=T)
P<-matrix(p,n.occasions-1,nrow=sum(marked),byrow=T)

#n.occasions<-dim(Phi)[2]+1
CH<-matrix(0,ncol=n.occasions,nrow=sum(marked))
#define a vector with marking occasion
mark.occ<-rep(1:length(marked),marked[1:length(marked)])
#fill in CH
for (i in 1:sum(marked))
      {
CH[i,mark.occ[i]]<-1
if (mark.occ[i]==n.occasions) next
   for(t in (mark.occ[i]+1):n.occasions)
         {
        #survive?
        sur<-rbinom(1,1,Phi[i,t-1])
         if(sur==0) break #move to next
         #recaptured?
         rp<-rbinom(1,1,P[i,t-1])
         if(rp==1) CH[i,t]<-1
           } #t
        } #i
return(CH)
}

#function to get occasion of marking for each animal
get.first<-function(x) min(which(x!=0))

###function to create capture history character strings (only need this for RMark runs)
pasty<-function(x) 
{
k<-ncol(x)
n<-nrow(x)
out<-array(dim=n)
for (i in 1:n)
{
out[i]<-paste(x[i,],collapse="")
}
return(out)
}

##################################



#FIRST MODEL: NO VARIATION OVER TIME IN PHI OR P

sink("cjs0.txt")
cat("
model {
#priors
mean.p~dunif(0,1)
mean.phi~dunif(0,1)

for (i in 1:nind)
 {
  for (t in f[i]:(n.occasions-1))
     {
      phi[i,t]<-mean.phi     #this allows us to expand in the future!
      p[i,t]<-mean.p
      }
    }
#likelihood
for (i in 1:nind)
 {
      z[i,f[i]]<-1   #state at first capture must be 1!
  for (t in (f[i]+1):n.occasions)
     {
        # state
	  mu1[i,t]<-phi[i,t-1]*z[i,t-1] 
         z[i,t]~dbern(mu1[i,t])
	    mu2[i,t]<-p[i,t-1]*z[i,t]
         y[i,t]~dbern(mu2[i,t])
         
         }
    }
}
",fill=TRUE)
sink()

#SECOND  MODEL: RANDOM TEMPORAL VARIATION OVER TIME IN PHI, CONSTANT P

sink("cjs_re.txt")
cat("
model {
#priors
mean.p~dunif(0,1)
mean.phi~dunif(0,1)
mu<-log(mean.phi/(1-mean.phi))
sigma~dunif(0,10)
tau<-pow(sigma,-2)
sigma2<-pow(sigma,2)

for (i in 1:nind)
 {
  for (t in f[i]:(n.occasions-1))

#PUTTING THINGS ON A LOGIT SCALE WILL MAKE THINGS MUCH EASIER LATER (E.G. FOR ADDING COVARIATES, OTHER HIERARCHICAL STRUCTURE)

     {logit(phi[i,t])<-mu+epsilon[t]

    
      p[i,t]<-mean.p
      }
    }

for (t in 1:(n.occasions-1))
 {
     epsilon[t]~dnorm(0,tau)
 }


#likelihood
for (i in 1:nind)
 {
      z[i,f[i]]<-1   #state at first capture must be 1!
  for (t in (f[i]+1):n.occasions)
     {
        # state
	  mu1[i,t]<-phi[i,t-1]*z[i,t-1] 
         z[i,t]~dbern(mu1[i,t])
	    mu2[i,t]<-p[i,t-1]*z[i,t]
         y[i,t]~dbern(mu2[i,t])
         
         }
    }



}
",fill=TRUE)
sink()

#THIRD  MODEL: FIXED TEMPORAL VARIATION OVER TIME IN PHI AND P

########################
sink("cjs_t.txt")
cat("
model {

#priors
  for (t in 1:(n.occasions-1))
     {
     phi[t]~dunif(0,1)
     p[t]~dunif(0,1)
    }



#likelihood
for (i in 1:nind)
 {
      z[i,f[i]]<-1   #state at first capture must be 1!
  for (t in (f[i]+1):n.occasions)
     {
        # state
	  mu1[i,t]<-phi[t-1]*z[i,t-1] 
         z[i,t]~dbern(mu1[i,t])
	    mu2[i,t]<-p[t-1]*z[i,t]
         y[i,t]~dbern(mu2[i,t])
         
         }
    }
}
",fill=TRUE)
sink()






#add to data known states (where we know z=1)
known.states.cjs<-function(ch){
state<-ch
for (i in 1:dim(ch)[1]){
n1<-min(which(ch[i,]==1))
n2<-max(which(ch[i,]==1))
state[i,n1:n2]<-1
state[i,n1]<-NA
 }
state[state==0]<-NA
return(state)
}

## initialize z states have to be NA where observed
cjs.init.z<-function(ch,f){
 for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
     n2<-max(which(ch[i,]==1))
     ch[i,f[i]:n2]<-NA
    }
 for (i in 1:dim(ch)[1])
  { ch[i,1:f[i]]<-NA
  }
 return(ch)
 }
#####SIMULATE DATA AS A RANDOM EFFECT
n.occas<-6
marked<-rep(50,n.occas-1)
mean.phi<-0.65
var.phi<-1

p.input<-rep(.4,n.occas-1)
logit.phi<-rnorm(n.occas-1,qlogis(mean.phi),var.phi^0.5)

#plogis is a built in inverse logit transformation
phi.input<-plogis(logit.phi)


CH<-simul.cjs(phi.input,p.input,marked)
f<-apply(CH,1,get.first)



#JAGS Models


jags.data <- list(y=CH,f=f,nind=dim(CH)[1],n.occasions=dim(CH)[2],z=known.states.cjs(CH))
####JAGS 
#homogenous p and phi)
params<-c("mean.phi","mean.p")
inits<-function(){list(z=cjs.init.z(CH,f),mean.phi=runif(1,0,1),mean.p=runif(1,0,1))}
m.0 <- jagsUI(data=jags.data, inits, parameters.to.save=params, model.file="cjs0.txt",n.thin = 5,n.chains = 3, n.burnin = 10000, n.iter =20000,parallel=TRUE)

summary(m.0)
# random effect model
params<-c("mean.phi","mean.p","sigma","phi")
inits<-function(){list(z=cjs.init.z(CH,f),mean.phi=runif(1,0,1),mean.p=runif(1,0,1),sigma=runif(1,0,10))}
m.re <- jagsUI(data=jags.data, inits, parameters.to.save=params, model.file="cjs_re.txt",n.thin = 5,n.chains = 3, n.burnin = 10000, n.iter =20000,parallel=TRUE)
summary(m.re)

#############
# ####fixed time effect model####
params<-c("phi","p")
inits<-function(){list(z=cjs.init.z(CH,f),phi=runif(n.occas-1,0,1),p=runif(n.occas-1,0,1),sigma=runif(1,0,10))}
m.t <- jagsUI(data=jags.data, inits, parameters.to.save=params, model.file="cjs_t.txt",n.thin = 5,n.chains = 3, n.burnin = 10000, n.iter =20000,parallel=TRUE)
summary(m.t)




dic.table<-data.frame(Model=c("Null","Random","Fixed"),DIC=c(m.0$DIC,m.re$DIC,m.t$DIC),pD=c(m.0$pD,m.re$pD,m.t$pD))

dic.table<-with(dic.table,dic.table[order(DIC),])

print("DIC comparison")
dic.table
#summaries output from models

print("Null Model")
summary(m.0)
print("Random Effect Model")
summary(m.re)
print("Fixed Effect Model")
summary(m.t)





#RMARK ANALYSIS for null and fixed time models
require(RMark)
CH.RMark<-pasty(CH)
rmark.data<-data.frame(ch=CH.RMark)
rmark.processed=process.data(rmark.data,model="CJS")
#default model is 'null'
m.0.mle<-mark(rmark.processed)
Phi.t=list(formula=~time)
p.t=list(formula=~time)

m.t.mle=mark(rmark.processed,model.parameters=list(Phi=Phi.t,p=p.t))
rm(list=ls())
cleanup(ask=FALSE)
