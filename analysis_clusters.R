# runs statistical analysis of population clusters in face space in guenon simulations
# Sandra Winters <sandra.winters@nyu.edu>

# setup ----------
library(data.table)
library(stringr)
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

nitt=2050000
burnin=50000
thin=500

# get data ----------
dat = read.csv('data_clusters.csv')

dat$simulation = str_replace(dat$simulation,'random_mating','RM')
dat$simulation = str_replace(dat$simulation,'mate_choice_PAM','PAM')
dat$simulation = str_replace(dat$simulation,'mate_choice','MC')

tmp = unlist(str_split(dat$simulation,'_',simplify=T)) 
dat$mating = tmp[,1]
dat$propAllo = as.numeric(str_replace(tmp[,2],'pAllo',''))
dat$hybridVia = as.numeric(str_replace(tmp[,3],'hybridVia',''))
dat$encFreq = as.numeric(str_replace(tmp[,4],'encFreq',''))
dat$npop = as.numeric(str_replace(tmp[,5],'populations',''))
rm(tmp)

dat$world = paste(dat$mating,dat$propAllo,dat$hybridVia,dat$encFreq,dat$npop,dat$iteration,sep='_')

dat$mating = factor(dat$mating)

dat$propAllo = dat$propAllo/100
dat$propSym = 1-dat$propAllo

dat$hybridVia = dat$hybridVia/100
dat$encFreq = dat$encFreq/100

dat.end = dat[dat$gen==20000,] # restrict to last generation

str(dat.end)

# plots ----------
hist(dat.end$propCorrect)

png('fig_clustering.png',width=22,height=7,units='cm',res=1200)
par(mfrow=c(1,3))
boxplot(propCorrect~mating,data=dat.end,ylab='Proportion of correctly classified faces',xlab='Mating pattern',outpch='.',boxwex=0.3,names=c('AMC','PAMC','NMC'))
# boxplot(propCorrect~propSym,data=dat.end,ylab='',xlab='Proportion of evolution in sympatry',outpch='.',boxwex=0.5)
boxplot(propCorrect~hybridVia,data=dat.end,ylab='',xlab='Hybrid fitness',outpch='.',boxwex=0.5)
# boxplot(propCorrect~encFreq,data=dat.end,ylab='',xlab='Conspecific encounter frequency',outpch='.',boxwex=0.5)
boxplot(propCorrect~npop,data=dat.end,ylab='',xlab='Number of co-evolving populations',outpch='.',boxwex=0.5)
dev.off()

par(mfrow=c(1,1))

# run models ----------
plotLoc = 'model_plots_clust/'
if (dir.exists(plotLoc)==0) {
  dir.create(plotLoc)
}

sink(paste0('model_results_clust.txt'))

for (i in 1:2) {
  if (i==1) {
    dat.end$mating = relevel(dat.end$mating,'RM')
    cat('\n---ORIGINAL MODEL---\n\n')
  } else {
    dat.end$mating = relevel(dat.end$mating,'MC')
    cat('\n---RELEVELED MODEL---\n\n')
  } 
  
  print(str(dat.end))
  
  m = MCMCglmm(c(nCorrect,nIncorrect)~mating+propSym+hybridVia+encFreq+npop,
               random=~world,
               family='multinomial2',
               prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=100))), #parameter expanded priors https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q3/004481.html
               data=dat.end,
               nitt=nitt,
               burnin=burnin,
               thin=thin,
               verbose=F)
  
  #diagnostics
  cat('\nMax autocorr in Sol:',max(autocorr(m$Sol)[2,,]),'\n')
  cat('\nMax autocorr in VCV:',max(autocorr(m$VCV)[2,,],na.rm=T),'\n')
  
  cat('\n')
  
  png(paste0(plotLoc,'m',i,'mean_sol1.png'))
  plot(unlist(m$Sol[,1:4]))
  dev.off()
  png(paste0(plotLoc,'m',i,'mean_sol2.png'))
  plot(unlist(m$Sol[,5:7]))
  dev.off()

  png(paste0(plotLoc,'m',i,'mean_vcv.png'))
  plot(m$VCV)
  dev.off()
  
  png(paste0(plotLoc,'m',i,'mean_autocorr.png'))
  plot.autocorr(m$Sol)
  dev.off()
  
  #stats
  m.coef <- t(coeftab(m)[,1])
  colnames(m.coef) <- rownames(coeftab(m))
  m.vcv <- var(m$Sol)
  
  cat('\n\nOverall ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=2:7)) #wald test for all fixed effects

  cat('\nMating ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=2:3)) #wald test for mating
  
  cat('\nProp sympatric ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=4)) #wald test for propSym
  
  cat('\nHybrid viability ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=5)) #wald test for hybridVia
  
  cat('\nEncounter frequency ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=6)) #wald test for encFreq
  
  cat('\nN populations ')
  print(wald.test(Sigma=m.vcv,b=as.vector(m.coef),Terms=7)) #wald test for npop
  
  print(summary(m))
  
  save.image(paste0('renv_clust_m',i,'.RData'))
}

sink()
file.show('model_results_clust.txt')

