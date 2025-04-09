#' Author: Emma Davies Smith, esmith at hsph.harvard.edu
#' Date: 23 OCT 2024

# Parameters -------------------------------------------------------------------

# Proportion of clusters allocated to treatment arm
r <- 0.5

# Total number of clusters
C <- c(10, 20, 30, 50)

# Average number of subjects per cluster
n <- c(10, 30)

# Cluster size distribution
n_type <- c('Equal', 'Unequal')

# Within-subject between-endpoint correlation
w12 <- c(-0.6, -0.5, -0.4)

# Intracluster correlations
p11 <- 0.10
p22 <- 0.05
p12 <- 0.025

# Endpoint win probabilities (null, small, moderate, large)
theta1 <- theta2 <- c(0.50, 0.56, 0.64, 0.71)

# Effect heterogeneity (difference in endpoint effects)
heterog <- c('None')

# Control parameters
pi0    <- 0.7    # Higher priority event probability
mu0    <- 50     # Lower priority mean
sigma  <- 5      # Lower priority SD (common to both arms)

# Number of simulation replicates
B <- 5000

# Scenarios --------------------------------------------------------------------

# Generates all scenarios based on input parameter values
scenarios <- expand.grid(
    B=B, r=r, C=C, n=n, n_type=n_type,
    w12=w12, p11=p11, p22=p22, p12=p12,
    theta1=theta1, theta2=theta2, heterog=heterog, 
    pi0=pi0, mu0=mu0, sigma=sigma
  ) %>% 
  # Effect heterogeneity
  mutate(theta2=case_when(heterog=='Large'    ~ theta1+0.105,
                          heterog=='Moderate' ~ theta1+0.07,
                          heterog=='Small'    ~ theta1+0.03,
                          heterog=='None'     ~ theta1)) 

# Attach treatment arm parameters req'd to obtain desired WinP
scenarios$pi1 <- mapply(FUN=function(x,y) getParam1(pi0=x, theta=y),
                        x=scenarios$pi0, y=scenarios$theta1)

scenarios$mu1 <- mapply(FUN=function(x,y,z) getParam2(mu0=x, sigma=y, theta=z),
                        x=scenarios$mu0, y=scenarios$sigma, z=scenarios$theta2)

# Attach true WinP
scenarios$theta <- mapply(FUN=function(a,b,c,d,e,f) 
                            getWinP(mu1=a, mu0=b, sigma=c, pi1=d, pi0=e, w12=f),
                          a=scenarios$mu1, b=scenarios$mu0, c=scenarios$sigma,
                          d=scenarios$pi1, e=scenarios$pi0, f=scenarios$w12)

# Arrange scenarios and add unique IDs
scenarios <- scenarios %>%
  arrange(r, n_type, heterog, theta, theta1, theta2, C, n, w12, pi0, mu0) %>%
  mutate(sid=row_number()) %>%
  dplyr::select(sid, r, n_type, heterog, theta, theta1, theta2, C, n, 
         w12, p11, p22, p12, pi0, pi1, mu0, mu1, sigma, B)

# Print total number of scenarios and preview
cat('Total scenarios:', nrow(scenarios), '\n')
head(scenarios)

# Output scenarios to dataframe for use in simulation
write.csv(scenarios, './out/scenarios.csv', row.names=F)
