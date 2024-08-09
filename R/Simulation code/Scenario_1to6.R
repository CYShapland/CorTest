#############################################################################################################################
# Testing Power of different estimators: Null Hypothesis that R=I and varying instruments
# Author: Chin Yang Shapland
# Last Updated: 28/09/22
############################################################################################################################

# Summary
# Checking the TIE from the correlation test from "psych" package,
# note that the covariance test is not here, as they are essentially testing the same thing,
# i.e. whether the off-diagonal is zero

### Set work directories ###
#wkdir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/identity_matrix/"
#setwd(wkdir)

### Load functions and packages ###
#fundir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/"
#source(paste(fundir, "Functions_CorTest_v6.R", sep=""))

### Set work directories ###
wkdir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/identity_matrix/"
setwd(wkdir)

### Load functions and packages ###
fundir<-"C:/Users/ew18103/University of Bristol/grp-Programme3-CYS - Documents/Bayesian_MissingData/R/CorrTests/"
source(paste(fundir, "Functions_CorTest_v6.R", sep=""))

library(MASS)
library(psych)
library("AER")

### Simulation ###
seed<-55
set.seed(seed)
nSim<-1000

#Parameters for relationship between G, X and Y
nSNP<-50
maf_min<-0.1
maf_max<-0.5
n<-8000
varXY<-0.1
varGX_range<-seq(0.05,0.45, by=0.2)

#Parameters for selection pressure
eta_0<-0.5
eta_z<-0
eta_c<-0
eta_y<-0
eta_x<-0.988

#Store results
res_nSNPdiff<-matrix(0, nSim, 6)
res_nSNPdiff_all<-list()

ptm <- proc.time()

for (k in 1:length(varGX_range)){
  #k<-1
  varGX<-varGX_range[k]

  for (i in 1:nSim){

    maf_dist<-runif(nSNP, maf_min, maf_max)

    sim_DataXY<-simData_contSNPs(maf_dist, (n*2), varXY, varGX)
    DataXY<-sim_DataXY$data

    #Standardising exposure and outcome
    DataXY$X<-(sim_DataXY$data[,"X"]-mean(sim_DataXY$data[,"X"]))/sd(sim_DataXY$data[,"X"])
    DataXY$Y<-(sim_DataXY$data[,"Y"]-mean(sim_DataXY$data[,"Y"]))/sd(sim_DataXY$data[,"Y"])

    #Randomly select individuals for missingness to ensure no overlap
    sample_MCAR<-sample(1:nrow(DataXY), n)
    sample_missX<-(1:nrow(DataXY))[-sample_MCAR]

    #Missing at random
    rand_select<-rbinom(n, 1, 0.6)
    DataXY_MCAR<-DataXY[sample_MCAR,]
    DataXY_MCAR$X<-ifelse(rand_select==1, DataXY_MCAR$X, NA)
    MCAR<-DataXY_MCAR[complete.cases(DataXY_MCAR), ]   # data randomly without missing X

    #Selection on X
    DataXY_missX<-DataXY[sample_missX,]
    lm_select<-eta_0+apply(eta_z*DataXY_missX[,paste("SNP",1:nSNP,sep="")],1, sum) + eta_x*DataXY_missX$X + eta_y*DataXY_missX$Y
    Prob_select<-exp(lm_select)/(1+exp(lm_select))
    DataXY_missX$X<-ifelse(Prob_select>0.6, DataXY_missX$X, NA)
    MissX<-DataXY_missX[complete.cases(DataXY_missX), ] # data with selection on X

    #Estimate correlation
    R_MCAR<-cor(MCAR[, paste("SNP",1:nSNP,sep="")])
    R_MissX<-cor(MissX[, paste("SNP",1:nSNP,sep="")])

    #the steiger test
    res_nSNPdiff[i,1]<-cortest(R_MCAR,NULL,n1=nrow(MCAR), n2=NULL)$prob
    res_nSNPdiff[i,2]<-cortest(R_MissX,NULL,n1=nrow(MissX), n2=NULL)$prob

    #the Jennrich test (n2=Inf avoids measurement error from identifity matrix)
    res_nSNPdiff[i,3]<-cortest_jennrich(R_MCAR,diag(nSNP),n1=nrow(MCAR), n2=Inf)$prob
    res_nSNPdiff[i,4]<-cortest_jennrich(R_MissX,diag(nSNP),n1=nrow(MissX), n2=Inf)$prob

    #an Bartlett test
    res_nSNPdiff[i,5]<-cortest.bartlett(R_MCAR,n=nrow(MCAR))$p.value
    res_nSNPdiff[i,6]<-cortest.bartlett(R_MissX,n=nrow(MissX))$p.value

  }

  res_nSNPdiff_all[[k]]<-res_nSNPdiff

}

RunTime<-proc.time() - ptm

etas<-c(Eta_0=eta_0, Eta_z=eta_z, Eta_c=eta_c, Eta_y=eta_y, Eta_x=eta_x)
names(res_nSNPdiff_all)<-varGX_range

saveResults(res_nSNPdiff_all, NA, "output/Res_CorTests_ident_VarGX", seed,nSim, nSNP, n, NA, varXY, varGX, etas, RunTime)
