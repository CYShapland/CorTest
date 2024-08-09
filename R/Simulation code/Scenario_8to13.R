#############################################################################################################################
# Two-sample testing
# Author: Chin Yang Shapland
# Last Updated: 16/12/22
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
nSim<-1000

#Parameters for relationship between G, X and Y
nSNP<-5
maf_min<-0.1
maf_max<-0.5
n_range<-seq(2000, 10000, by=2000)
varXY<-0.1
varGX<-0.9

#Parameters for selection pressure
eta_0<-0.5
eta_z<-0
eta_c<-0
eta_y<-0
eta_x<-0.988

#Store results
res_nSNPdiff_all<-list()
SimCheck_all<-list()

ptm <- proc.time()

for (k in 1:length(n_range)){
  #k<-1
  n<-n_range[k]

  SimCheck<-matrix(0, nSim, 2)
  res_nSNPdiff<-matrix(0, nSim, 6)

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
    lm_select<-eta_0+apply(eta_z*DataXY_missX[,paste("SNP",1:nSNP,sep="")],1, sum) + eta_x*DataXY_missX$X + eta_y*DataXY_missX$Y
    Prob_select<-exp(lm_select)/(1+exp(lm_select))
    DataXY_missX$X<-sapply(1:n, function(x) ifelse(Prob_select[x]>=runif(1, 0, 1), DataXY_missX$X[x], NA))
    MissX_1<-DataXY_missX[complete.cases(DataXY_missX), ] # data with selection on X
    MissX_2<-DataXY_missX[!complete.cases(DataXY_missX), ] # data with missing X

    #Checking simulation does give mean of 60% selection with sd of 0.2
    SimCheck[i,1]<-mean(Prob_select)
    SimCheck[i,2]<-sd(Prob_select)

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
