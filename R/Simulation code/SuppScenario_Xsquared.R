#############################################################################################################################
# Testing power of different estimators: Null Hypothesis that R_1=R_2
# Author: Chin Yang Shapland
# Last Updated: 28/09/22
############################################################################################################################

# Summary
# Checking the power from the correlation test from "psych" package, and also BoXM test

### Set work directories ###
wkdir<-"/mnt/storage/home/ew18103/CorTest/TwoSample/"
setwd(wkdir)

### Load functions and packages ###
source("/mnt/storage/home/ew18103/CorTest/Functions_CorTest_v6.R")

library(MASS)
library(psych)
library("AER")
library(heplots)

### Simulation ###
seed<-55
set.seed(seed)
nSim<-1000

#Parameters for relationship between G, X and Y
nSNP<-50
maf_min<-0.1
maf_max<-0.5
n_range<-seq(2000, 10000, by=2000)
varXY<-0.1
varGX<-0.45

#Parameters for selection pressure
eta_0<-0.5
eta_z<-0
eta_c<-0
eta_y<-0
eta_x<-0.988

#Store results
SimCheck_all<-list()
res_oneSample_all<-list()

ptm <- proc.time()

for (k in 1:length(varGX_range)){
  #k<-1
  varGX<-varGX_range[k]

  #Store results
  SimCheck<-matrix(0, nSim, 4)
  res_oneSample<-matrix(0, nSim, 12)

  for (i in 1:nSim){

    maf_dist<-runif(nSNP, maf_min, maf_max)

    sim_DataXY<-simData(maf_dist, (n*2), varXY, varGX)
    DataXY<-sim_DataXY$data

    #Standardising exposure and outcome
    DataXY$X<-(sim_DataXY$data[,"X"]-mean(sim_DataXY$data[,"X"]))/sd(sim_DataXY$data[,"X"])

    #Randomly select individuals for missingness to ensure no overlap
    sample_MCAR<-sample(1:nrow(DataXY), n)
    sample_missX<-(1:nrow(DataXY))[-sample_MCAR]

    #Missing at random
    DataXY_MCAR<-DataXY[sample_MCAR,]
    rand_select<-rbinom(n, 1, 0.6)
    DataXY_MCAR$X<-ifelse(rand_select==1, DataXY_MCAR$X, NA)
    MCAR_1<-DataXY_MCAR[complete.cases(DataXY_MCAR), ]   # data randomly without missing X
    MCAR_2<-DataXY_MCAR[!complete.cases(DataXY_MCAR), ]  # data with missing X

    #Selection on X
    DataXY_missX<-DataXY[sample_missX,]
    lm_select<-eta_0 + eta_x*(DataXY_missX$X)^2
    Prob_select<-exp(lm_select)/(1+exp(lm_select))
    DataXY_missX$X<-sapply(1:n, function(x) ifelse(Prob_select[x]>=runif(1, 0, 1), DataXY_missX$X[x], NA))
    MissX_1<-DataXY_missX[complete.cases(DataXY_missX), ] # data with selection on X
    MissX_2<-DataXY_missX[!complete.cases(DataXY_missX), ] # data with missing X

    #Checking simulation does give mean of 60% selection with sd of 0.2
    SimCheck[i,1]<-mean(Prob_select)
    SimCheck[i,2]<-sd(Prob_select)
    SimCheck[i,3]<-nrow(MissX_1)
    SimCheck[i,4]<-nrow(MCAR_1)

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

    #the Jennrich test (n2=Inf avoids measurement error from identity matrix)
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

etas<-c(Eta_0=eta_0, Eta_z=eta_z, Eta_c=eta_c, Eta_y=eta_y, Eta_x=eta_x)
names(res_oneSample_all)<-varGX_range

saveResults(res_oneSample_all, SimCheck_all, "output/Res_CovCorTests_ident_N", seed,nSim, nSNP, n, NA, varXY, varGX, etas, RunTime)
