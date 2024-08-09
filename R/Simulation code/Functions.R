#########################################################################################################
# Function required simulate different scenarios
#########################################################################################################

my_packages<-c("gmm", "ivmodel", "MASS")
lapply(my_packages, require, character.only = TRUE)

###################################################
# Simulation of X with N(0,1) and continuous SNP
###################################################

simData_normX <- function(n_SNP, n, XYVar, ZXVar) {

  ### Genotype ###
  z <- matrix(0,n,n_SNP)
  for(i in 1: n_SNP ) {
    z[,i] <- rnorm(n,0,1)
  }

  ### X has variance of 1 ###
  u <- rnorm(n,0,1)
  rx<-rnorm(n,0,1)

  #Their coefficient to explain 2% of X is
  alpha<-sqrt(0.5*(1-ZXVar))

  ZXVar_zi<-ZXVar/n_SNP
  beta_norm<-sqrt(ZXVar_zi)

  mx<-0

  X_norm<-mx + rowSums(beta_norm*z) + alpha*u + alpha*rx

  geno <- data.frame(z,X_norm)
  colnames(geno) <-c(paste("z",1:n_SNP,sep=""), "X")

  return(geno)
}

######################
# Simulation of SNPs
######################

simData <- function(MAF, n, XYVar, GXVar) {

  n_SNP <- length(MAF)

  ### Genotype ###
  z <- matrix(0,n,n_SNP)
  for(i in 1: n_SNP ) {
    z[,i] <- rnorm(n,0,1)
  }

  cor.SNP<-cor(z)

  SNPs<-matrix(0,n,n_SNP)

  for (r in 1:n_SNP) {
    p1 <- qnorm((1-MAF[r])^2)
    p2 <- qnorm(1-MAF[r]^2)
    SNPs[,r]<-ifelse(z[,r] >= p1,1,0)
    SNPs[,r]<-ifelse(z[,r] >= p2,2,SNPs[,r])
  }

  ### X has variance of 1 ###
  u <- rnorm(n,0,1)
  rx<-rnorm(n,0,1)

  #Their coefficient to explain 2% of X is
  alpha<-sqrt(0.5*(1-GXVar))

  GXVar_SNPi<-GXVar/n_SNP
  beta<-sqrt(GXVar_SNPi/(2*MAF*(1-MAF)))

  mx<-0

  X<-mx + rowSums(t(beta*t(SNPs))) + alpha*u + alpha*rx

  ### Y has variance of 1 ###

  #confounder and random coefficient to explain 6% of Y is
  ry<-rnorm(n,0,1)

  delta<-sqrt(0.5*(1-XYVar))

  phi<-sqrt(XYVar)

  my<-0

  Y<-my + phi*X + delta*u + delta*ry

  geno <- data.frame(SNPs,X,Y)
  colnames(geno)[1:(n_SNP+2)] <-c(paste("SNP",1:n_SNP,sep=""), "X", "Y")

  returnlist<-list(data=geno, beta_xy=phi,cor.snp=cor.SNP, maf=MAF)
  return(returnlist)
}
