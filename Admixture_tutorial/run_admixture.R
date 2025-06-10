###############   Admixture ########################

setwd('/Users/osipova/Documents/LabDocs/TA_human_PopGen/Admixture_tutorial/')

ADMIXTURE <- '/Users/osipova/Documents/LabDocs/TA_human_PopGen/software/admixture_macosx-1.3.0/admixture '




## run Admixture for k = 2
k <- 2
system(paste0(ADMIXTURE,
              ' hapmap_r23_subset.bed ',
              k))

## Read output
runs <- list()
runs[[k]]<-read.table(paste0("hapmap_r23_subset.", k, ".Q"))
barplot(t(as.matrix(runs[[k]])), col=rainbow(k), ylab="Ancestry", border=NA)



## run Admixture for k = 1,2,3,4,5
for (k in 1:5){
  system(paste0(ADMIXTURE,
                '--cv ',
                'hapmap_r23_subset.bed ',
                k,
                ' > ', k, '.log'))
              }


## Load Admixture output
runs <- list()
for (k in 1:5){
  runs[[k]]<-read.table(paste0("hapmap_r23_subset.", k, ".Q"))
}

## Plot runs for all k
par(mfrow=c(4,1))
for (k in 2:5){
  barplot(t(as.matrix(runs[[k]])), col=rainbow(k), ylab="Ancestry", border=NA)
}

