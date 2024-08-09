#############################################################################################################################
# Two-sample Testing with continuous variables
# Author: Chin Yang Shapland
# Last Updated: 11/01/24
############################################################################################################################

### Load functions and packages ###
source("Functions.R")

library(MASS)
library(psych)
library("AER")
library(heplots)

### Simulation ###
seed<-55
set.seed(seed)
nSim<-10

#Parameters for relationship between G, X and Y
nSNP<-50
n<-100000
varXY<-0.1
varZX<-0.05

#Parameters for selection pressure
eta_0<-0.5
eta_z<-0
eta_c<-0
eta_y<-0
eta_x<-0.988

#Checks that Prob_Select has mean of 0.6 and sd of 0.2
SimCheck_all<-list()

#Store results
res_nSNPdiff_all<-list()

for (k in 1:length(n_range)){
  #k<-1
  n<-n_range[k]

  #store results
  SimCheck<-matrix(0, nSim, 4)
  res_nSNPdiff<-matrix(0, nSim, 6)

  ptm <- proc.time()

  for (i in 1:nSim){
    print(i)

    DataXY<-simData_normX(nSNP, n*2, varXY, varZX)

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
    res_nSNPdiff[i,5]<-boxM(DataXY_MCAR[,paste("z",1:nSNP,sep="")], DataXY_MCAR[,"Miss"])$p.value
    res_nSNPdiff[i,6]<-boxM(DataXY_missX[,paste("z",1:nSNP,sep="")], DataXY_missX[,"Miss"])$p.value

  }

  RunTime<-proc.time() - ptm

  res_nSNPdiff_all[[k]]<-res_nSNPdiff
  SimCheck_all[[k]]<-SimCheck
}
