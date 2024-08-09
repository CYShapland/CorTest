#########################################################################################################
# Function required for testing different correlation comparison tests
# Last Updated: 05/01/24
# Summary:
# Extracted from "SettingUp_manyIVs_060521.R" and "Functions_CorrTest_Zheng2019.R"
#########################################################################################################

# Update:
# (1) Add simulation to continuous SNPs
# (2) function for doing 2SLS and
# (3) extracting results from 2SLS outputs
# (4) Update iv models packages
# (5) MAF function
# (6) Apostolos corrected Jennrich
# (7) Adding function from Luepsen
# (8) Added correlated SNPs
# (9) Speeding up function from Luepsen using cor()
# (10) Simulate X with N(0,1) and continuous SNP
# (11) Saving both OneSample and TwoSample Tests

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

#################################
# Simulation of continuous SNPs
#################################

simData_contSNPs <- function(MAF, n, XYVar, GXVar) {

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
  #GXVar<-0.6
  alpha<-sqrt(0.5*(1-GXVar))

  beta<-sqrt((GXVar/n_SNP)/(2*MAF*(1-MAF)))
  mx<-0

  X<-mx + rowSums(t(beta*t(SNPs))) + alpha*u + alpha*rx

  ### Y has variance of 1 ###

  #confounder and random coefficient to explain 6% of Y is
  ry<-rnorm(n,0,1)

  delta<-sqrt(0.5*(1-XYVar))

  phi<-sqrt(XYVar)

  my<-0

  Y<-my + phi*X + delta*u + delta*ry

  geno <- data.frame(z,SNPs,X,Y)
  colnames(geno) <-c(paste("z",1:n_SNP,sep=""), paste("SNP",1:n_SNP,sep=""), "X", "Y")

  returnlist<-list(data=geno, beta_xy=phi,cor.snp=cor.SNP, maf=MAF)
  return(returnlist)
}

#####################
# Correlation Tests #
#####################

# Apostolos modified Jennrich test
# See email 23/09/2022 "Manuscript: Sensitivity analysis to detect selection on exposure using genetic correlation"
# Make sure you use n2=Inf if you want to do one-sample to avoid measurement error

cortest_jennrich <- function (R1, R2, n1 = NULL, n2 = NULL) {

  p <- dim(R1)[2]
  if (dim(R1)[1] != p) {
    n1 <- dim(R1)[1]
    R1 <- cor(R1, use = "pairwise")
    warning("R1 matrix was not square, correlations found")
  }
  if (dim(R2)[1] != dim(R2)[2]) {
    n2 <- dim(R2)[1]
    R2 <- cor(R2, use = "pairwise")
    warning("R2 matrix was not square, correlations found")
  }
  if (!is.matrix(R1))
    R1 <- as.matrix(R1)
  if (!is.matrix(R2))
    R2 <- as.matrix(R2)
  if (dim(R1)[2] != dim(R2)[2])
    stop("correlation matrices M and S must be of the same size!")
  if (is.null(n2))
    n2 <- n1

  if (!is.null(n1) & !is.null(n2)) {
    if (n2 == Inf) {
      c <- n1
      R <- R2
    } else {
      c <- n1 * n2/(n1 + n2)
      R <- (n1 * R1 + n2 * R2)/(n1 + n2)
    }
  } else {
    c <- 1
    R <- (n1 * R1 + n2 * R2)/(n1 + n2)
  }

  R.inv <- solve(R)
  S <- diag(p) + R * R.inv
  S.inv <- solve(S)
  R.diff <- R1 - R2
  Z <- sqrt(c) * R.inv %*% R.diff
  chi2 <- tr(Z %*% t(Z))/2 - t(diag(Z)) %*% S.inv %*% diag(Z)
  chi2 <- chi2[1, 1]
  df <- p * (p - 1)/2
  results <- list(chi2 = chi2, prob = pchisq(chi2, df, lower.tail = FALSE), df=df)
  return(results)

}

#################
# Saving output #
#################

### saveResults() ###

saveResults<-function(result_OneSample, result_TwoSample, Sim_checks, file, seed,nsim, n_SNP, n, rmax, varXY, varGX, etas, RunTime){
  choices<-c(seed,nsim,n_SNP,n, rmax, varXY, varGX, etas, RunTime)
  names(choices)<-c("seed","nsim","n_SNP", "sample_size", "rmax","R2_XY", "R2_GX", "Eta_0", "Eta_z", "Eta_c",
                    "Eta_y", "Eta_x", "user", "system", "elapsed")

  save(result_OneSample, result_TwoSample, Sim_checks,choices, file=paste(file, ".Rdata",sep=""))
}

###################
# MAF Calculation #
###################

### maf_cal() ###

maf_cal<-function(x){
  X<-as.vector(table(x))
  X <- X/sum(X)
  p <- (X[1] + 0.5 * X[2])/sum(X)
  y <- pmin(p, 1 - p)
  print(y)
}

