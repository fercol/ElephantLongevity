# ================================== CODE METADATA =============================
# AUTHOR: Fernando Colchero
# DATE CREATED: 2023-12-14
# DATE MODIFIED: 
# DESCRIPTION: Functions for reproducibility code for Sach et al. (2024)
# NOTES: Uses packages:
#           - paramDemo (parametric and non-parametric demographic functions)
#           - BaSTA2.0 (Bayesian survival trajectory analysis)
#           - snowfall (parallel computing)
#           - mvtnorm (multivariate truncated normal)
# ================================ CODE START ==================================
# ==================== #
# ==== LIBRARIES: ====
# ==================== #
# --------------------------- #
# ---- Install packages: ----
# --------------------------- #
# Your installed packages:
instPacks <- installed.packages()[, 1]

# RColorBrewer:
if (!"RColorBrewer" %in% instPacks) {
  install.packages("RColorBrewer")
}

# Devtools to install GitHub packages:
if (!"devtools" %in% instPacks) {
  install.packages("devtools")
}

# mvtnorm:
if (!"mvtnorm" %in% instPacks) {
  install.packages("mvtnorm")
}

# paramDemo:
if (!"paramDemo" %in% instPacks) {
  devtools::install_git("https://github.com/fercol/paramDemo", subdir = "pkg/")
}

# install BaSTA2.0:
if (!"BaSTA2.0" %in% instPacks) {
  devtools::install_git("https://github.com/fercol/BaSTA2.0", subdir = "pkg/")
}

# snowfall for parallel computing:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}

# ------------------------ #
# ---- Load packages: ----
# ------------------------ #
library(mvtnorm)
library(paramDemo)
library(BaSTA2.0)
library(RColorBrewer)

# =============================== #
# ==== ADDITIONAL FUNCTIONS: ==== 
# =============================== #
# Kulback-Leibler discrepancy:
CalcKLc <- function(mean1, mean2, sd1, sd2, low = -Inf) {
  pMat <- cbind(p1 = c(mean = mean1, sd = sd1),
                p2 = c(mean = mean2, sd = sd2))
  parRan <- range(sapply(1:2, function(pp) {
    qtnorm(c(0.001, 0.999), pMat["mean", pp], pMat["sd", pp], lower = low)
  }))
  parVec <- seq(parRan[1], parRan[2], length = 100)
  dp <- parVec[2] - parVec[1]
  parDens <- sapply(1:2, function(pp) 
    dtnorm(parVec, pMat["mean", pp], pMat["sd", pp], lower = low))
  p1dens <- parDens[, 1]
  p2dens <- parDens[, 2]
  idp <- which(p1dens > 0 & p2dens > 0)
  kld <- sum(p2dens[idp] * log(p2dens[idp] / p1dens[idp]) * dp)
  qKlc <- (1 + (1 - exp(-2 * kld)^(1 / 2))) / 2
  outList <- c(KL = kld, qKL = qKlc)
  return(outList)
}

# Truncated normal functions:
dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# Calculate range:
CalcRan <- function(cxy, mxy, alpha = 0.001) {
  ran <- sapply(1:2, function(ixy) {
    qnorm(p = c(alpha / 2, 1 - alpha / 2), mean = mxy[ixy], 
          sd = sqrt(cxy[ixy, ixy]))
  })
  colnames(ran) <- c("x", "y")
  return(ran)  
}

# Calculate bivariate density:
CalcBivDens <- function(cxy, mxy, ranxy, n = 100) {
  xyv <- sapply(1:2, function(ixy) {
    seq(ranxy[1, ixy], ranxy[2, ixy], length = n)
  })
  colnames(xyv) <- c("x", "y")
  xycom <- cbind(x = rep(xyv[, "x"], nrow(xyv)), y = rep(xyv[, "y"], 
                                                         each = nrow(xyv)))
  
  dxy <- xyv[2, ] - xyv[1, ]
  dd <- dmvnorm(xycom, mean = mxy, sigma = cxy)
  res <- list(dens = dd, dxy = dxy)
  return(res)
}

# Calculate bivariate Kullback-Liebler:
CalcBivKL <- function(m1, m2, cv1, cv2, n = 100) {
  ran1 <- CalcRan(cv1, m1)
  ran2 <- CalcRan(cv2, m2)
  
  ranxy <- ran1
  for (ii in 1:2) {
    ranxy[1, ii] <- min(ran1[1, ii], ran2[1, ii])
    ranxy[2, ii] <- max(ran1[2, ii], ran2[2, ii])
  }
  
  dens1 <- CalcBivDens(cxy = cv1, mxy = m1, ranxy = ranxy, n = n)
  dens2 <- CalcBivDens(cxy = cv2, mxy = m2, ranxy = ranxy, n = n)
  idn0 <- which(dens1$dens > 0 & dens2$dens > 0)
  
  kl <- sum(dens2$dens[idn0] * log(dens2$dens[idn0] / dens1$dens[idn0]) * 
              prod(dens1$dxy))
  
  qKl <- ((1 + (1 - exp(-2 * kl)^(1 / 2))) / 2 - 0.5) * 2
  return(c(KL = kl, qKL = qKl))
}
