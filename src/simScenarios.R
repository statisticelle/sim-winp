#' Simulate scenarios
#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024
#'
#' @param scenarios dataframe containing all scenarios to simulate
#' @param out name (string) of csv file to output (to out subdir)
#'
#' @return Replicate summary statistics by scenario
simScenarios <- function(scenarios, out, seed){
  
  # Set random seed
  set.seed(seed)
  
  # Initialize dataframe to store results for each simulation scenario
  scenResults <- data.frame(
    sid=NA,      # scenario id
    fail=NA,     # number of times lmm failed to converge (across reps)
    est=NA,      # mean point estimate (across reps)
    se=NA,       # mean standard error
    icc=NA,      # mean intracluster correlation
    identL=NA,   # identity ci: mean lower limit
    identU=NA,   # identity ci: mean upper limit
    identC=NA,   # identity ci: empirical coverage probability (ecp)
    identR=NA,   # identity ci: empirical rejection rate (err)
    identEL=NA,  # identity ci: empirical left tail error rate
    identER=NA,  # identity ci: empirical right tail error rate
    logitL=NA,   # logit ci: mean lower limit
    logitU=NA,   # logit ci: mean upper limit
    logitC=NA,   # logit ci: empirical coverage probability (ecp)
    logitR=NA,   # logit ci: empirical rejection rate (err)
    logitEL=NA,  # logit ci: empirical left tail error rate
    logitER=NA   # logit ci: empirical right tail error rate
  )
  
  # Total number of scenarios to investigate
  n_scen <- max(scenarios$sid)
  
  for (s in 1:n_scen){
    cat('SCENARIO:', s, '\n')
    
    # Extract relevant scenario
    scenario <- scenarios[s,]
    theta    <- scenario$theta # true parameter value
    
    # Get corresponding target correlation matrix R
    R  <- getCorrelation(n=scenario$n,
                         w12=scenario$w12,
                         p11=scenario$p11,
                         p22=scenario$p22,
                         p12=scenario$p12)
    
    cat('Generating data for all replicates...\n')
    # Generate control responses on original scale for all replicates
    C0 <- (1-scenario$r)*scenario$C
    X0 <- genResponse(C=C0, 
                      n=scenario$n,
                      n_type=scenario$n_type,
                      B=scenario$B,
                      R=R,
                      Pi=scenario$pi0,
                      Mu=scenario$mu0,
                      Sigma=scenario$sigma) 
    X0$Arm <- 0 
    
    # Generate treatment responses on original scale for all replicates
    C1 <- (scenario$r)*scenario$C
    X1 <- genResponse(C=C1, 
                      n=scenario$n,
                      n_type=scenario$n_type,
                      B=scenario$B,
                      R=R,
                      Pi=scenario$pi1,
                      Mu=scenario$mu1,
                      Sigma=scenario$sigma)
    X1$Arm <- 1 
    
    # Combine and label all responses, assign unique cluster labels for modelling
    X <- rbind(X0, X1) %>% 
      arrange( B, Arm, C, pid) %>%
      group_by(B, Arm, C) %>%
      mutate(cid=cur_group_id()) %>%
      ungroup() %>%
      mutate(pid=row_number()) 
    
    # Initialize dataframe to store results for each scenario replicate
    repResults <- data.frame(
      b=NA,        # replicate number
      est=NA,      # point estimate
      se=NA,       # standard error
      icc=NA,      # intracluster correlation
      identL=NA,   # identity ci: lower limit
      identU=NA,   # identity ci: upper limit
      identC=NA,   # identity ci: coverage? yes=1, no=0
      identR=NA,   # identity ci: reject null? yes=1, no=0
      identEL=NA,  # identity ci: left tail error? yes=1, no=0
      identER=NA,  # identity ci: right tail error? yes=1, no=0
      logitL=NA,   # logit ci: lower limit
      logitU=NA,   # logit ci: upper limit
      logitC=NA,   # logit ci: coverage? yes=1, no=0
      logitR=NA,   # logit ci: reject null? yes=1, no=0
      logitEL=NA,  # logit ci: left tail error? yes=1, no=0
      logitER=NA   # logit ci: right tail error? yes=1, no=0
    )

    
    cat('Fitting and assessing model for each replicate...\n')
    fail <- 0 # Track convergence failures
    for (b in 1:scenario$B){
      if (b %% 1000 == 0) cat('REPLICATE:', b, '\n')
      
      # Dataframe for replicate b
      Xb <- X[X$B==b,]
      Yb <- getFractions(Xb)
      
      # Fit linear mixed model to win fractions
      # Using optim for more robust optimization
      mod <- try(
        lme(YP~Arm, random=~1|cid, data=Yb, 
            control=lmeControl(opt='optim'))
      )
      
      if (!inherits(mod, "try-error")){ # If no convergence failure...
        # Fixed effect estimates and variances
        mod_fest <- fixef(mod)
        mod_fvar <- vcov(mod)
        mod_fdf  <- mod$fixDF$X
        
        # Random effect variances
        mod_rvar <- matrix(as.numeric(VarCorr(mod)), ncol=2)
        
        # Construct win probability estimates
        est <- (mod_fest[2]+1)/2    # winp = (beta1+1)/2
        se  <- sqrt(mod_fvar[2,2])  # se(winp) = se(beta1)
        icc <- mod_rvar[1,1]/(mod_rvar[1,1] + mod_rvar[2,1])
        
        alpha <- 0.05
        t <- qt(1-alpha/2, mod_fdf[2])
        
        # Identity confidence interval
        ident_ci <- est + c(-1,1)*t*se
        identC   <- checkCoverage(theta, ident_ci)
        identR   <- checkReject(ident_ci)
        identEL  <- checkLeft(theta, ident_ci)
        identER  <- checkRight(theta, ident_ci)
        
        # Logit confidence interval
        logit_lu <- logit(est) + c(-1,1)*t*se/(est*(1-est))
        logit_ci <- expit(logit_lu)
        logitC   <- checkCoverage(theta, logit_ci)
        logitR   <- checkReject(logit_ci)
        logitEL  <- checkLeft(theta, logit_ci)
        logitER  <- checkRight(theta, logit_ci)
        
        # Store results for replicate b
        result <- c(b, est, se, icc, 
                    ident_ci[1], ident_ci[2], identC, identR, identEL, identER,
                    logit_ci[1], logit_ci[2], logitC, logitR, logitEL, logitER)
        repResults[b,] <- result
      } else { # If convergence failure...
        # Increment fail count by 1
        fail <- fail + 1
      }
      
    }
    
    # Average across all replicates within scenario and store summary result
    repSummary <- colMeans(repResults, na.rm=T)[-1]
    scenResults[s,] <- c('sid'=s, 'fail'=fail/scenario$B, repSummary)
    path <- paste0('./out/', out) 
    write.csv(scenResults, path, row.names=F)
    
  }
  return(scenResults)
}