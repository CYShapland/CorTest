#############################################################################################################################
# Testing power of different estimators: Null Hypothesis that R_1=R_2
# Author: Chin Yang Shapland
# Last Updated: 29/02/24
############################################################################################################################

# Summary
# Had to change threshold approach to make sure the 60% sample selection is still true.
# used truncated normal to get mean of 0.6 with sd of 0.2

### Set work directories ###
wkdir<-"/mnt/storage/home/ew18103/CorTest/TwoSample/"
setwd(wkdir)

### Load functions and packages ###
source("/mnt/storage/home/ew18103/CorTest/Functions_CorTest_v6.R")

library(MASS)
library(psych)
library("AER")
library(heplots)
library(tidyverse)
library(truncnorm)

### Simulation ###
seed<-55
set.seed(seed)
nSim<-1000

#Parameters for relationship between G, X and Y
nSNP<-50
maf_min<-0.1
maf_max<-0.5
n<-10000
varXY<-0.1
varGX<-0.45

#Parameters for selection pressure
thresX_range<-seq(0.15,0.5, by=0.1)

#Store results
SimCheck_all<-list()
res_nSNPdiff_all<-list()

ptm <- proc.time()

for (k in 1:length(thresX_range)){
  #k<-1
  thresX<-thresX_range[k]

  #Store results
  SimCheck<-matrix(0, nSim, 2)
  res_nSNPdiff<-matrix(0, nSim, 6)

  for (i in 1:nSim){

    maf_dist<-runif(nSNP, maf_min, maf_max)

    sim_DataXY<-simData(maf_dist, (n*5), varXY, varGX)
    DataXY<-sim_DataXY$data

    #Standardising exposure and outcome
    DataXY$X<-(sim_DataXY$data[,"X"]-mean(sim_DataXY$data[,"X"]))/sd(sim_DataXY$data[,"X"])

    #Randomly select individuals for missingness to ensure no overlap
    sample_MCAR<-sample(1:nrow(DataXY), n)

    #Missing at random
    DataXY_MCAR<-DataXY[sample_MCAR,]
    rand_select<-rbinom(n, 1, 0.6)
    DataXY_MCAR$X<-ifelse(rand_select==1, DataXY_MCAR$X, NA)
    MCAR_1<-DataXY_MCAR[complete.cases(DataXY_MCAR), ]   # data randomly without missing X
    MCAR_2<-DataXY_MCAR[!complete.cases(DataXY_MCAR), ]  # data with missing X

    #Selection on X
    repeat{
      SelectProp<-rtruncnorm(1, a=0, b=1, mean = 0.6, sd = 0.2)
      DataXY_missX <- DataXY %>%
        mutate(X = if_else(X > quantile(X, thresX), X, NA))
      if ((n*(1-SelectProp))<sum(is.na(DataXY_missX$X))) break
    }
    MissX_1<-DataXY_missX[sample(which(!is.na(DataXY_missX$X)),round((n*SelectProp))), ] # data with selection on X
    MissX_2<-DataXY_missX[sample(which(is.na(DataXY_missX$X)),round(n*(1-SelectProp))), ] # data with missing X

    #Checking simulation does give mean of 60% selection with sd of 0.2
    SimCheck[i,1]<-nrow(MissX_1)/n
    SimCheck[i,2]<-nrow(MCAR_1)/n

    #Estimate correlation
    MCAR_R1<-cor(MCAR_1[,paste("SNP",1:nSNP,sep="")])
    MCAR_R2<-cor(MCAR_2[,paste("SNP",1:nSNP,sep="")])

    MissX_R1<-cor(MissX_1[,paste("SNP",1:nSNP,sep="")])
    MissX_R2<-cor(MissX_2[,paste("SNP",1:nSNP,sep="")])

    #the steiger test
    res_nSNPdiff[i,1]<-cortest(MCAR_R1,  MCAR_R2, n1=nrow(MCAR_1), n2=nrow(MCAR_2))$prob
    res_nSNPdiff[i,2]<-cortest(MissX_R1, MissX_R2, n1=nrow(MissX_1), n2=nrow(MissX_2))$prob

    #the Jennrich test
    res_nSNPdiff[i,3]<-cortest_jennrich(MCAR_R1,  MCAR_R2, n1=nrow(MCAR_1), n2=nrow(MCAR_2))$prob
    res_nSNPdiff[i,4]<-cortest_jennrich(MissX_R1, MissX_R2, n1=nrow(MissX_1), n2=nrow(MissX_2))$prob

    #add missing indicator for BoxM
    DataXY_MCAR["Miss"]<-ifelse(is.na(DataXY_MCAR$X), 1, 2)
    DataXY_missX["Miss"]<-ifelse(is.na(DataXY_missX$X), 1, 2)

    #BoxM test
    res_nSNPdiff[i,5]<-boxM(DataXY_MCAR[,paste("SNP",1:nSNP,sep="")], DataXY_MCAR[,"Miss"])$p.value
    res_nSNPdiff[i,6]<-boxM(DataXY_missX[,paste("SNP",1:nSNP,sep="")], DataXY_missX[,"Miss"])$p.value


  }
  res_nSNPdiff_all[[k]]<-res_nSNPdiff
  SimCheck_all[[k]]<-SimCheck
}

RunTime<-proc.time() - ptm

etas<-rep(NA,5)
names(res_nSNPdiff_all)<-thresX_range

saveResults(res_nSNPdiff_all, SimCheck_all, "output/Res_CovCorTests_TwoSample_ThresX_v2", seed,nSim, nSNP, n, NA, varXY, varGX, etas, RunTime)
