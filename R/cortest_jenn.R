cortest_jennrich <- function (R1, R2, n1 = NULL, n2 = NULL) {

  R1<-MCAR_R1
  R2<-MCAR_R2
  n1=as.numeric(nrow(MCAR_1))
  n2=as.numeric(nrow(MCAR_2))

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

  MCAR_R1,  MCAR_R2
  n1=as.numeric(nrow(MCAR_1))
  n2=as.numeric(nrow(MCAR_2))

  N<-n1+n2
  rbar<-((n1 * MCAR_R1) + (n2 * MCAR_R2))/N
  rbar.inv <- solve(rbar)
  Z1<-sqrt(n1) * rbar.inv %*% (MCAR_R1-rbar)
  Z2<-sqrt(n2) * rbar.inv %*% (MCAR_R2-rbar)

  S <- diag(nSNP) + (rbar %*% rbar.inv)
  S.inv <- solve(S)

  chi2_1<-tr(Z1 %*% t(Z1))/2 - t(diag(Z1)) %*% S.inv %*% diag(Z1)
  chi2_2<-tr(Z2 %*% t(Z2))/2 - t(diag(Z2)) %*% S.inv %*% diag(Z2)

  chi2_1+chi2_2

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
