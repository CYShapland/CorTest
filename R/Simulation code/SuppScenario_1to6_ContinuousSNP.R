#############################################################################################################################
# One-sample Testing with continuous variables
# Author: Chin Yang Shapland
# Last Updated: 28/09/22
############################################################################################################################

### Load functions and packages ###
source("Functions.R")

library(MASS)
library(psych)
library("AER")

### Simulation ###
seed<-55
set.seed(seed)
nSim<-1000

#Parameters for relationship between G, X and Y
nSNP<-50
n_range<-seq(2000, 10000, by=2000)
varXY<-0.1
varZX<-0.45

#Parameters for selection pressure
eta_0<-0.5
eta_z<-0
eta_c<-0
eta_y<-0
eta_x<-0.988

#Store results
res_oneSample_all<-list()
SimCheck_all<-list()

for (k in 1:length(n_range)){
  #k<-1
  n<-n_range[k]

  #store result
  res_oneSample<-matrix(0, nSim, 12)
  SimCheck<-matrix(0, nSim, 4)

  for (i in 1:nSim){

    DataXY<-simData_normX(nSNP, n*2, varXY, varZX)

    #Standardising exposure and outcome
    DataXY$X<-(DataXY[,"X"]-mean(DataXY[,"X"]))/sd(DataXY[,"X"])

    #Randomly select individuals for missingness to ensure no overlap
    sample_MCAR<-sample(1:nrow(DataXY), n)
    sample_missX<-(1:nrow(DataXY))[-sample_MCAR]

    #Missing at random
    rand_select<-rbinom(n, 1, 0.6)
    DataXY_MCAR<-DataXY[sample_MCAR,]
    DataXY_MCAR$X<-ifelse(rand_select==1, DataXY_MCAR$X, NA)
    MCAR_1<-DataXY_MCAR[complete.cases(DataXY_MCAR), ]   # data randomly without missing X
    MCAR_2<-DataXY_MCAR[!complete.cases(DataXY_MCAR), ]  # data with missing X

    #Selection on X
    DataXY_missX<-DataXY[sample_missX,]
    lm_select<-eta_0 + eta_x*DataXY_missX$X
    Prob_select<-exp(lm_select)/(1+exp(lm_select))
    DataXY_missX$X<-sapply(1:n, function(x) ifelse(Prob_select[x]>=runif(1, 0, 1), DataXY_missX$X[x], NA))
    MissX_1<-DataXY_missX[complete.cases(DataXY_missX), ] # data with selection on X
    MissX_2<-DataXY_missX[!complete.cases(DataXY_missX), ] # data with missing X

    #Checking simulation does give mean of 60% selection with sd of 0.2
    SimCheck[i,1]<-mean(Prob_select)
    SimCheck[i,2]<-sd(Prob_select)
    SimCheck[i,3]<-nrow(MissX_1)/n
    SimCheck[i,4]<-nrow(MCAR_1)/n

    #Estimate correlation
    MCAR_R1<-cor(MCAR_1[,paste("z",1:nSNP,sep="")])
    MCAR_R2<-cor(MCAR_2[,paste("z",1:nSNP,sep="")])

    MissX_R1<-cor(MissX_1[,paste("z",1:nSNP,sep="")])
    MissX_R2<-cor(MissX_2[,paste("z",1:nSNP,sep="")])

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
