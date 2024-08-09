#############################################################################################################################
# Testing power of different estimators: Null Hypothesis that R_1=R_2
# Author: Chin Yang Shapland
# Last Updated: 28/09/22
############################################################################################################################

# Summary
# Checking the power from the correlation test from "psych" package, and also BoXM test

### Set work directories ###
wkdir<-"/mnt/storage/home/ew18103/CorTest/TwoSample/" #bc4
#wkdir<-"/user/home/ew18103/CovCorTest/OneSample"  #bluepebble
setwd(wkdir)

### Load functions and packages ###
source("/mnt/storage/home/ew18103/CorTest/Functions_CorTest_v6.R") #bc4
#source("/user/home/ew18103/CovCorTest/Functions_CorTest_v6.R")  #bluepebble

### Set work directories ###
#wkdir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/identity_matrix/ProbSelect/"
#setwd(wkdir)

### Load functions and packages ###
#fundir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/"
#source(paste(fundir, "Functions_CorTest_v6.R", sep=""))

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
res_oneSample_all<-list()

ptm <- proc.time()

for (k in 1:length(thresX_range)){
  #k<-1
  thresX<-thresX_range[k]

  #Store results
  SimCheck<-matrix(0, nSim, 2)
  res_oneSample<-matrix(0, nSim, 12)

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
    table(is.na(DataXY_MCAR$X))

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

    ### One sample test ###
    #the steiger test
    res_oneSample[i,1]<-cortest(MCAR_R1,NULL,n1=nrow(MCAR_1), n2=NULL)$prob
    res_oneSample[i,2]<-cortest(MCAR_R2,NULL,n1=nrow(MCAR_2), n2=NULL)$prob
    res_oneSample[i,3]<-cortest(MissX_R1,NULL,n1=nrow(MissX_1), n2=NULL)$prob
    res_oneSample[i,4]<-cortest(MissX_R2,NULL,n1=nrow(MissX_2), n2=NULL)$prob

    #the Jennrich test (n2=Inf avoids measurement error from identifity matrix)
    res_oneSample[i,5]<-cortest_jennrich(MCAR_R1,diag(nSNP),n1=nrow(MCAR_1), n2=Inf)$prob
    res_oneSample[i,6]<-cortest_jennrich(MCAR_R2,diag(nSNP),n1=nrow(MCAR_2), n2=Inf)$prob
    res_oneSample[i,7]<-cortest_jennrich(MissX_R1,diag(nSNP),n1=nrow(MissX_1), n2=Inf)$prob
    res_oneSample[i,8]<-cortest_jennrich(MissX_R2,diag(nSNP),n1=nrow(MissX_2), n2=Inf)$prob

    #an Bartlett test
    res_oneSample[i,9]<-cortest.bartlett(MCAR_R1,n=nrow(MCAR_1))$p.value
    res_oneSample[i,10]<-cortest.bartlett(MCAR_R2,n=nrow(MCAR_2))$p.value
    res_oneSample[i,11]<-cortest.bartlett(MissX_R1,n=nrow(MissX_1))$p.value
    res_oneSample[i,12]<-cortest.bartlett(MissX_R2,n=nrow(MissX_2))$p.value

  }
  res_oneSample_all[[k]]<-res_oneSample
  SimCheck_all[[k]]<-SimCheck
}

RunTime<-proc.time() - ptm

etas<-rep(NA,5)
names(res_oneSample_all)<-thresX_range

saveResults(res_oneSample_all, SimCheck_all, "output/Res_CovCorTests_identMatrix_ThresX_v2", seed,nSim, nSNP, n, NA, varXY, varGX, etas, RunTime)
