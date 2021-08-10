# Investigate:
# Easy ("Good"): Look at variance of Wp^p (wasserstein p distance to the p) and Wp when ?? = ?? as N increases
# Loop, store distance, compute estimator of sample variance (1/n ? sum((xi - x??)^2))

# var(c(10 9 8 7 6 5 4 3 2 1))
# pattern: want plot variance against N, definintely for Wp^p goes down faster than 1/n, interesting: find out the exponent.
# it's n^-c c> 1 but how big is this and how does it depend on dimension
# How to see what size of c is:
# Imagine I know a-priori it's log(var) = logc * c log(n)
# log log plot of variance against sample size then slope = c < -1
# Then take p^th root then calculate variance and then loglog plot of variance against sample size
# Pattern should be the same for whatever metric and exponent p > 1
# Which exact metric to take and which exponent to take is pretty wide
# Start with euclidean norm exponent 2, then integer exponents
# Sample size up to 10,000
# Replicates: A couple

# Wasserstein takes a metric d and then (integral( d(x - y) ^p)^1/p

# d <- function(x, y, q = 2) {(sum(x - y)^q)^(1/q)} # = Lq-norm which then induces a metric but yeah

# preliminary investigation: CLT papers on Optimal transport ?? if I would like (only, otherwise do not, i repeat do not unless you would like)
# Came out february, good papers -> you must enjoy them. contains refernece to previous work which is pretty good
# https://arxiv.org/pdf/2102.06379.pdf

set.seed(1)

library(tidyverse)
library(glue)
library(Rcpp)
library(RcppEigen)
sourceCpp("src/network_simplex_fast.cpp") # Compile and load the solver. The header files are also necessary.
sourceCpp("src/costMatrix.cpp")
sourceCpp("src/squared_cost_eval.cpp")

n <- 100  # Sample size
d <- 1   # Dimension
p <- 2    # Wasserstein-p norm
q <- 2    # Wasserstein-p norm
K <- 100 # Number of replicates
R <- 1000 # Number of bootstrap replicates
ns <- c(100, 200, 300)
# ns <- c(100, 300, 1000, 3000, 10000)
ds <- c(1, 10, 100)   # Dimension
var_wp_emp <- rep(NA, length(ns))
var_wp_emp_p <- rep(NA, length(ns))
wp_boot_sd <- rep(NA, length(ns))
wp_boot_p_sd <- rep(NA, length(ns))

j <- 0
for (n in ns){
    wp_emp_p <- rep(NA, K)
    j <- j+1
    for (k in 1:K){
        # 1.) Generate samples
        muhat <- matrix(rnorm(n * d), nrow = d)
        nuhat <- matrix(rnorm(n * d), nrow = d)
        
        # 2.) Compute the cost matrix 
        cost_mat <- costMatrix(muhat, nuhat, q, p) # d(x, y) = L-q metric, Wasserstein-p norm
        
        # 3.) Compute the distance from the cost matrix
        # K samples, variance of wasserstein-p distance over 100 samples
        wp_emp_p[k] <- (SolveAssignmentNetworkflow(n, cost_mat))
    }
    
    wp_emp <- wp_emp_p^1/p
    
    # 4.) Compute the variance
    var_wp_emp[j] <- var(wp_emp)
    var_wp_emp_p[j] <- var(wp_emp_p)
    
    # 5.) Calculate bootstrap distribution
    
    wp_boot <- rep(NA, R)
    wp_boot_p <- rep(NA, R)
    for (i in 1:R){
        sample <- sample(wp_emp, K, replace = T)
        wp_boot[i] <- var(sample)   
        wp_boot_p[i] <- var(sample^p)   
    }
    
    # 6.) Calculate confidence intervals of original estimate
    wp_boot_sd[j] <- sd(wp_boot)
    wp_boot_p_sd[j] <- sd(wp_boot_p)
}

wp_emp_lci <- var_wp_emp - 2*wp_boot_sd
wp_emp_uci <- var_wp_emp + 2*wp_boot_sd

# wp_emp_lci <- var_wp_emp_p - 2*wp_boot_p_sd
# wp_emp_uci <- var_wp_emp_p + 2*wp_boot_p_sd

# Plot variance against sample size
# Sample size = ns, variance = var_wp_emp
# Also plot wp_emp_lcu and uci in blue or smth idk

coefs <- lm(log10(var_wp_emp) ~ log10(ns))$coefficients

qplot(ns, var_wp_emp) +
    theme_bw() +
    scale_x_log10() +
    scale_y_log10() +
    geom_ribbon(aes(x = ns, ymin = wp_emp_lci, ymax = wp_emp_uci), alpha = 0.3, fill = "#008ccc") +
    geom_abline(intercept = coefs[1], slope = coefs[2], linetype = "dashed", size = 0.7)

# Same plot with wp_p ^

# Loop over a few ds, same plots, separate d = 1 case and d = 10 case
# test with p = 

# increase d, the slope of the line should get worse, variance may decrease slower
# 100% Correct CONJECTURE: in dim-1 whatever the norm is, and whatever p is, slope should be -2
# extra fun: prove the conjecture above
# Continue with this I guess: nevermind, there's a bit of maths there
