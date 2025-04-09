#' Helper functions
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

#' Logit transformation
logit <- function(x){
  log(x/(1-x))
}

#' Expit transformation (logit back-transformation)
expit <- function(x){
  exp(x)/(1+exp(x))
}

#' Check true value of parameter falls within confidence limits
checkCoverage <- function(theta, ci){
  ifelse(theta >= ci[1] & theta <= ci[2], 1, 0)
}

#' Check null value (0.5) falls outside confidence limits
checkReject <- function(ci){
  ifelse(ci[1] <= 0.5 & ci[2] >= 0.5, 0, 1)
}

# Check if truth falls to left of lower limit
checkLeft <- function(theta, ci){
  ifelse(theta <= ci[1], 1, 0)
}

# Check if truth falls to right of upper limit
checkRight <- function(theta, ci){
  ifelse(theta >= ci[2], 1, 0)
}

# Calculate mean of continuous endpoint conditional on binary event
cMean <- function(x, mu, sigma, pi, w12){
  val <- mu + sigma / sqrt(pi*(1-pi)) * w12 * (x-pi)
  return(val)
}

# Calculate variance of continuous endpoint conditional on binary event
cVar <- function(sigma, w12){
  val <- sigma^2 * (1-w12^2)
  return(val)
}

# Calculate true win probability given parameters for both arms
getWinP <- function(mu1, mu0, sigma, pi1, pi0, w12){
  
  # mean within treatment conditional on event/non-event
  mu11 <- cMean(x=1, mu=mu1, sigma=sigma, pi=pi1, w12=w12)
  mu10 <- cMean(x=0, mu=mu1, sigma=sigma, pi=pi1, w12=w12)
  
  # mean within control conditional on event/non-event
  mu01 <- cMean(x=1, mu=mu0, sigma=sigma, pi=pi0, w12=w12)
  mu00 <- cMean(x=0, mu=mu0, sigma=sigma, pi=pi0, w12=w12)
  
  # conditional variance (same for all strata)
  # sigma11=sigma10=sigma01=sigma00=cvar
  cvar <- cVar(sigma=sigma, w12=w12)
  
  # calculate standardized mean differences among event/non-event
  smd1 <- (mu11-mu01)/sqrt(2*cvar)
  smd0 <- (mu10-mu00)/sqrt(2*cvar)
  
  # calculate win probabilities among event/non-event
  theta21 <- pnorm(smd1) 
  theta20 <- pnorm(smd0)
  
  # calculate win probability for binary outcome (highest priority)
  theta1 <- 0.5 + 0.5*(pi0-pi1)
  
  # calculate true value of win probability parameter
  theta <- theta1 + (theta21-0.5)*pi0*pi1 + (theta20-0.5)*(1-pi0)*(1-pi1)
  
  return(theta)
}

#' Construct target block exhangeable correlation matrix
#' @details
#' For more details on block exchangeable correlation matrix see:
#' Wang et al. (CCT 2021) A flexible sample size solution for longitudinal and
#' crossover cluster randomized trials with continuous outcomes.
#' 
#' @param n   Number of subjects per cluster
#' @param w12 Within-subject between-endpoint correlation
#' @param p11 Intracluster correlation of 1st endpoint
#' @param p22 Intracluster correlation of 2nd endpoint
#' @param p12 Within-cluster between-subject between-endpoint correlation
#'
#' @return Target correlation matrix R
#'
#' @examples
#' getCorrelation(n=3, w12=0.3, p11=0.10, p22=0.05, p12=0.025)
getCorrelation <- function(n, w12, p11, p22, p12){
  
  # Endpoint correlation
  Omega <- matrix(c(1, w12, w12, 1), ncol=2)
  
  # Intracluster correlations
  Phi <- matrix(c(p11, p12, p12, p22), ncol=2)
  
  # Unit matrices
  I <- diag(1, nrow=n, ncol=n)
  J <- matrix(1, nrow=n, ncol=n)
  
  # Desired block exchangeable correlation matrix
  R <- I %x% (Omega-Phi) + J %x% Phi
  return(R)
}

#' Get treatment arm parameters for first endpoint given control param and WinP
#' @details
#' Binary endpoint with pi0 > pi1. Win probability is function of risk 
#' difference.
#' 
#' @param pi0 Event probability in control arm (fixed)
#' @param theta Desired win probability 
#'
#' @return Event probability in treatment arm
#'
#' @examples
#' getParam1(pi0=0.3, theta=0.5)  # return should be equal to pi0
#' getParam1(pi0=0.7, theta=0.64) # return should be equal to 0.42
getParam1 <- function(pi0, theta){
  pi1 <- pi0 + 1 - 2*theta 
  return(pi1)
}

#' Get treatment arm parameters for second endpoint given control param and WinP
#' @details
#' Normal endpoint with mu1 > mu0. Win probability is function of standardized
#' mean difference.
#' 
#' @param mu0 Mean in control arm (fixed)
#' @param sigma SD shared by control and treatment arms (fixed)
#' @param theta Desired win probability 
#'
#' @return Mean in treatment arm.
#'
#' @examples
#' getParam2(mu0=50, sigma=5, theta=0.5)  # return should be equal to mu0
#' getParam2(mu0=50, sigma=5, theta=0.64) # return should be equal to ~52.5
getParam2 <- function(mu0, sigma, theta){
  mu1 <- mu0 + sqrt(2*sigma^2) * qnorm(theta)
  return(mu1)
}

getFractions <- function(X){
  # X is data frame containing bivariate responses (X1, X2) for Arm=0 and Arm=1
  
  # Create all N0xN1 pairs of control-treatment participants
  P1 <- unique(X$pid[X$Arm==1])     # Distinct treatment pids
  P0 <- unique(X$pid[X$Arm==0])     # Distinct control pids
  
  pairs <- expand.grid('P1'=P1, 'P0'=P0) %>%                # Create all pairs
    left_join(X, by=c('P1'='pid')) %>%                      # Attach trt resp
    dplyr::select(P1, P0, 'X11'=X1, 'X12'=X2) %>%         
    left_join(X, by=c('P0'='pid')) %>%                      # Attach ctrl resp
    dplyr::select(P1, P0, X11, X12, 'X01'=X1, 'X02'=X2)
  
  # Construct prioritized Heaviside scores
  scores <- pairs %>%
    mutate(s1 = (X11<X01)-(X11>X01), # Treat non-event vs ctrl event => win
           s2 = (X12>X02)-(X12<X02), # Treat greater than ctrl => win
           sP = s1 + s2*(s1==0),     # Prioritized sign score
           hP = (1+sP)/2)            # Prioritized Heaviside score
  
  # Create win fractions by arm and combine into one dataset
  Y0 <- scores %>% 
    group_by(P0) %>%
    summarize(YP=mean(1-hP)) %>%
    dplyr::select('pid'=P0, YP)
  
  Y1 <- scores %>%
    group_by(P1) %>%
    summarize(YP=mean(hP)) %>%
    dplyr::select('pid'=P1, YP)
  
  Y <- rbind(Y0, Y1) %>%
    left_join(X, by='pid')
  
  return(Y)
}

