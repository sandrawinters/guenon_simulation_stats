# compiles data from guenon simulations for statistical analyses
# Sandra Winters <sandra.winters@nyu.edu>

# setup ----------
library(data.table)
library(stringr)
library(tidyr)

resDir = getwd() #SET TO LOCATION OF 'simulation_data' FOLDER

# get simulation data (compiles face_data files) ----------
dat = data.frame(mating=character(),
                 propAllo=numeric(),
                 hybridVia=numeric(),
                 encFreq=numeric(),
                 npop=numeric(),
                 itt=numeric(),
                 gen=numeric(),
                 pop=numeric(),
                 meanDistCon=numeric(),
                 sdDistCon=numeric(),
                 meanDistHetero=numeric(),
                 sdDistHetero=numeric(),
                 meanQuality=numeric(),
                 sdQuality=numeric(),
                 meanBias=numeric(),
                 sdBias=numeric(),
                 stringsAsFactors=F)

for (m in c('mate_choice','random_mating','mate_choice_PAM')) { 
  for (a in c(10,50)) {
    for (h in c(0,2,5,10,50,90)) {
      for (e in c(25,50,75)) {
        for (p in 2:6) {
          simName = paste0(m,'_',
                           as.character(a),'pAllo_',
                           as.character(h),'hybridVia_',
                           as.character(e),'encFreq_',
                           as.character(p),'populations')
          cat(simName,'\n')
          d = read.csv(paste0(resDir, 
                              '/simulation_data/',
                              simName,
                              '_face_data.csv'))
          d$mating = m
          d$propAllo = a
          d$hybridVia = h
          d$encFreq = e
          d$npop = p

          dat = rbind(dat,d)
        }
      }
    }
  }
}

# save.image('tmp.RData')

#restrict to final generation
dat = droplevels(dat[dat$gen==20000,])

#format data
setcolorder(dat,neworder='npop')
setcolorder(dat,neworder='encFreq')
setcolorder(dat,neworder='hybridVia')
setcolorder(dat,neworder='propAllo')
setcolorder(dat,neworder='mating')

dat$mating = str_replace(dat$mating,'random_mating','RM')
dat$mating = str_replace(dat$mating,'mate_choice_PAM','PAM')
dat$mating = str_replace(dat$mating,'mate_choice','MC')

dat$world = paste(dat$mating,dat$propAllo,dat$hybridVia,dat$encFreq,dat$npop,dat$itt,sep='_')

# get face distance data & format ----------
dat.dist = read.csv('data_distance.csv')
dat.dist = dat.dist[,c('sim','itt','pop','stat','spp','gen20000')]

colnames(dat.dist)[colnames(dat.dist)=='gen20000'] = 'distance'

dat.dist$sim = str_replace(dat.dist$sim,'random_mating','RM')
dat.dist$sim = str_replace(dat.dist$sim,'mate_choice_PAM','PAM')
dat.dist$sim = str_replace(dat.dist$sim,'mate_choice','MC')

dat.dist$stat[dat.dist$stat==1] = 'mean'
dat.dist$stat[dat.dist$stat==2] = 'std'

dat.dist$spp[dat.dist$spp==1] = 'conspp'
dat.dist$spp[dat.dist$spp==2] = 'heterospp'

tmp = unlist(str_split(dat.dist$sim,'_',simplify=T))
dat.dist$mating = tmp[,1]
dat.dist$propAllo = as.numeric(str_replace(tmp[,2],'pAllo',''))
dat.dist$hybridVia = as.numeric(str_replace(tmp[,3],'hybridVia',''))
dat.dist$encFreq = as.numeric(str_replace(tmp[,4],'encFreq',''))
dat.dist$npop = as.numeric(str_replace(tmp[,5],'populations',''))
rm(tmp)

dat.dist$world = paste(dat.dist$mating,dat.dist$propAllo,dat.dist$hybridVia,dat.dist$encFreq,dat.dist$npop,dat.dist$itt,sep='_')

# calculate mean distance between population average faces, mean distance within population faces ----------
#restrict distance data to mean distance to conspecific faces
dat.dist = dat.dist[dat.dist$stat=='mean' & dat.dist$spp=='conspp',]

#distance function
euclidean <- function(a, b) sqrt(sum((a - b)^2))

#calculate mean distances 
dat$meanDistBT = NaN; #between population average faces (for each population, calculate mean pairwise distance between it and all others it co-evolved with)
dat$meanDistWI = NaN; #within faces within the same population
for (i in 1:dim(dat)[1]) {
  #get mean faces for all other populations in the simulation
  d = dat[dat$world==dat$world[i] & dat$pop!=dat$pop[i],grep('faceMean',colnames(dat))]
  
  #calculate average distance to other population mean faces
  sumDist = 0;
  for (j in 1:dim(d)[1]) {
    sumDist = sumDist + euclidean(d[j,],dat[i,grep('faceMean',colnames(dat))])
  }
  dat$meanDistBT[i] = sumDist/dim(d)[1]
  
  #get mean distance between dyads within population
  dat$meanDistWI[i] = dat.dist$distance[dat.dist$world==dat$world[i] & dat.dist$pop==dat$pop[i]]
}

#save
write.csv(dat,'simulation_data.csv',row.names=F)
