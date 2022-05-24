# runs statistical analysis of evolved female mating biases in guenon simulations
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

nitt=2050000
burnin=50000
thin=500

# get data ----------
dat = read.csv('data_variation_end.csv')

dat = dat[,c('mating','propAllo','hybridVia','encFreq','npop','itt','pop','biasMean','meanDistBT','meanDistWI')]

dat$mating = str_replace(dat$mating,'random_mating','RM')
dat$mating = str_replace(dat$mating,'mate_choice_PAM','PAM')
dat$mating = str_replace(dat$mating,'mate_choice','MC')
dat$mating = factor(dat$mating)

dat$world = paste(dat$mating,dat$propAllo,dat$hybridVia,dat$encFreq,dat$npop,dat$itt,sep='_')

dat$propAllo = dat$propAllo/100
dat$propSym = 1-dat$propAllo

dat$hybridVia = dat$hybridVia/100
dat$encFreq = dat$encFreq/100

# plots ----------
hist(dat$biasMean)

par(mfrow=c(1,5))
boxplot(biasMean~mating,data=dat,xlab='Mating',ylab='Likelihood of engaging in mate choice',boxwex=0.5)
boxplot(biasMean~propSym,data=dat,xlab='Proportion of evolution in sympatry',ylab='',boxwex=0.5)
boxplot(biasMean~hybridVia,data=dat,xlab='Hybrid fitness',ylab='',boxwex=0.5)
boxplot(biasMean~encFreq,data=dat,xlab='Conspecific encounter frequency',ylab='',boxwex=0.5)
boxplot(biasMean~npop,data=dat,xlab='Number of coevolving populations',ylab='',boxwex=0.5)
par(mfrow=c(1,1))

# run model ----------
plotLoc = 'model_plots_pref/'
if (dir.exists(plotLoc)==0) {
  dir.create(plotLoc)
}

sink(paste0('model_results_pref.txt'))

print(str(dat))

#model
m = MCMCglmm(biasMean~mating+propSym+hybridVia+encFreq+npop,
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

png(paste0(plotLoc,'mean_sol1.png'))
plot(unlist(m$Sol[,1:4]))
dev.off()
png(paste0(plotLoc,'mean_sol2.png'))
plot(unlist(m$Sol[,5:7]))
dev.off()

png(paste0(plotLoc,'mean_vcv.png'))
plot(m$VCV)
dev.off()

png(paste0(plotLoc,'mean_autocorr.png'))
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

save.image(paste0('renv_pref.RData'))

sink()
file.show('model_results_pref.txt')


# compare preferences to population distances ----------
#run models
plotLoc = 'model_plots_pref_dist/'
if (dir.exists(plotLoc)==0) {
  dir.create(plotLoc)
}

sink(paste0('model_results_pref_dist.txt'))

print(str(dat))

for (i in 1:2) {
  if (i==1) {
    modelID = 'distBT'
    dat$dist = dat$meanDistBT
    cat('\n===BIAS ~ MEAN DISTANCE BETWEN POPULATIONS===\n\n')
  } else {
    modelID = 'distWI'
    dat$dist = dat$meanDistWI
    cat('\n---BIAS ~ MEAN DISTANCE WITHIN POPULATIONS===\n\n')
  }

  #model
  m = MCMCglmm(biasMean~dist,
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
  
  png(paste0(plotLoc,modelID,'_mean_sol.png'))
  plot(unlist(m$Sol))
  dev.off()
  
  png(paste0(plotLoc,modelID,'_mean_vcv.png'))
  plot(m$VCV)
  dev.off()
  
  png(paste0(plotLoc,modelID,'_mean_autocorr.png'))
  plot.autocorr(m$Sol)
  dev.off()
  
  #stats
  cat('\n----------\nModel summary\n----------\n ')
  print(summary(m))
  
  save.image(paste0('renv_pref_dist_',modelID,'.RData'))
}

sink()
file.show('model_results_pref_dist.txt')

