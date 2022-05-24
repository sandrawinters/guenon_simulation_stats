# runs statistical analysis of evolved population distances in face space in guenon simulations
# Sandra Winters <sandra.winters@nyu.edu>

# setup ----------
library(tidyr)
library(MCMCglmm)
library(aod)
library(stringr)
library(coefplot2)
#coefplot2 has problems with R-Forge build -- to install: 
#(1) install packages reshape & lme4
#(2) install coefplot2 from Bolker's website: 
#install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source")

source('coeftab.R')

plot.autocorr <- function(x) {
  n <- dim(x)[2]
  par(mfrow=c(ceiling(n/2),2), mar=c(3,2,3,0))
  for (i in 1:n) {
    acf(x[,i], lag.max=100, main=colnames(x)[i])
    grid()
  }
}

nitt=160000
burnin=10000
thin=50

# get data ----------
dat = read.csv('data_variation_end.csv',stringsAsFactors=F)

dat$mating = factor(dat$mating)

dat$propAllo = dat$propAllo/100
dat$propSym = 1-dat$propAllo

dat$hybridVia = dat$hybridVia/100;
dat$encFreq = dat$encFreq/100;

# plots ----------
hist(dat$meanDistBT)

par(mfrow=c(1,5))
boxplot(meanDistBT~mating,data=dat,ylab='Mean distance between populations',xlab='Mating pattern',outpch='.',boxwex=0.3,names=c('AMC','PAMC','NMC'))
boxplot(meanDistBT~propSym,data=dat,ylab='',xlab='Proportion of evolution in sympatry',outpch='.',boxwex=0.5)
boxplot(meanDistBT~hybridVia,data=dat,ylab='',xlab='Hybrid fitness',outpch='.',boxwex=0.5)
boxplot(meanDistBT~encFreq,data=dat,ylab='',xlab='Conspecific encounter frequency',outpch='.',boxwex=0.5)
boxplot(meanDistBT~npop,data=dat,ylab='',xlab='Number of co-evolving populations',outpch='.',boxwex=0.5)
par(mfrow=c(1,1))

hist(dat$meanDistWI)

par(mfrow=c(1,5))
boxplot(meanDistWI~mating,data=dat,ylab='Mean distance within populations',xlab='Mating pattern',outpch='.',boxwex=0.3,names=c('AMC','PAMC','NMC'))
boxplot(meanDistWI~propSym,data=dat,ylab='',xlab='Proportion of evolution in sympatry',outpch='.',boxwex=0.5)
boxplot(meanDistWI~hybridVia,data=dat,ylab='',xlab='Hybrid fitness',outpch='.',boxwex=0.5)
boxplot(meanDistWI~encFreq,data=dat,ylab='',xlab='Conspecific encounter frequency',outpch='.',boxwex=0.5)
boxplot(meanDistWI~npop,data=dat,ylab='',xlab='Number of co-evolving populations',outpch='.',boxwex=0.5)
par(mfrow=c(1,1))

# run models ----------
plotLoc = 'model_plots/'
if (dir.exists(plotLoc)==0) {
  dir.create(plotLoc)
}

sink(paste0('model_results_dist.txt'))

for (i in 1:4) {
  if (i==1) {
    modelID = 'distBT1'
    dat$dist = dat$meanDistBT
    dat$mating = relevel(dat$mating,'RM')
    cat('\n===MEAN DISTANCE BETWEN POPULATIONS===\n\n')
  } else if (i==2) {
    modelID = 'distBT2'
    dat$dist = dat$meanDistBT
    dat$mating = relevel(dat$mating,'MC')
    cat('\n===MEAN DISTANCE BETWEEN POPULATIONS - RELEVELED===\n\n')
  } else if (i==3) {
    modelID = 'distWI1'
    dat$dist = dat$meanDistWI
    dat$mating = relevel(dat$mating,'RM')
    cat('\n===MEAN DISTANCE WITHIN POPULATIONS===\n\n')
  } else {
    modelID = 'distWI2'
    dat$dist = dat$meanDistWI
    dat$mating = relevel(dat$mating,'MC')
    cat('\n===MEAN DISTANCE WITHIN POPULATIONS - RELEVELED===\n\n')
  } 
  
  print(str(dat))
  
  #model
  m = MCMCglmm(dist~mating+propSym+hybridVia+encFreq+npop+mating*hybridVia,
               random=~world,
               family='gaussian',
               prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002))),
               data=dat,
               nitt=nitt,
               burnin=burnin,
               thin=thin,
               verbose=F)
  
  #diagnostics
  cat('\nMax autocorr in Sol:',max(autocorr(m$Sol)[2,,]),'\n')
  cat('\nMax autocorr in VCV:',max(autocorr(m$VCV)[2,,],na.rm=T),'\n')
  
  cat('\n')
  
  png(paste0(plotLoc,modelID,'_sol1.png'))
  plot(unlist(m$Sol[,1:4]))
  dev.off()
  png(paste0(plotLoc,modelID,'_sol2.png'))
  plot(unlist(m$Sol[,5:7]))
  dev.off()

  png(paste0(plotLoc,modelID,'_vcv.png'))
  plot(m$VCV)
  dev.off()
  
  png(paste0(plotLoc,modelID,'_autocorr.png'))
  plot.autocorr(m$Sol)
  dev.off()
  
  #stats
  m.coef <- t(coeftab(m)[,1])
  colnames(m.coef) <- rownames(coeftab(m))
  m.vcv <- var(m$Sol)
  
  cat('\n\n----------\nOverall ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=2:7)) #wald test for all fixed effects

  cat('\n----------\nMating ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=2:3)) #wald test for mating
  
  cat('\n----------\nProp sympatric ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=4)) #wald test for propSym
  
  cat('\n----------\nHybrid viability ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=5)) #wald test for hybridVia
  
  cat('\n----------\nEncounter frequency ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=6)) #wald test for encFreq
  
  cat('\n----------\nN populations ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=7)) #wald test for npop
  
  cat('\n----------\nModel summary\n----------\n ')
  print(summary(m))
  print(HPDinterval(mcmc(m$Sol)))
  
  save.image(paste0('renv_',modelID,'.RData'))
}
  
sink()
file.show('model_results_dist.txt')

