set.seed(1)
library(tidyverse)
library(Rcpp)
library(RcppEigen)
sourceCpp("src/network_simplex_fast.cpp") # Compile and load the solver. The header files are also necessary.
sourceCpp("src/squared_cost_eval.cpp")

n <- 100  # Sample size
d <- 1   # Dimension
p <- 2    # Wasserstein-p norm
K <- 1000 # Number of replicates
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
        
        # wp_emp_p sort, diff, mean, it's somewhere bro just like find it
        
        # 2.) Compute the cost matrix 
        cost_mat <- EvaluateSquaredCost(muhat, nuhat, 4L) # d(x, y) = L-q metric, Wasserstein-p norm
        
        # 3.) Compute the distance from the cost matrix
        # K samples, variance of wasserstein-p distance over 100 samples
        wp_emp_p[k] <- SolveAssignmentNetworkflow(n, cost_mat)
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

coefs[2]

# Same plot with wp_p ^

# Loop over a few ds, same plots, separate d = 1 case and d = 10 case
# test with p = 

# increase d, the slope of the line should get worse, variance may decrease slower
# 100% Correct CONJECTURE: in dim-1 whatever the norm is, and whatever p is, slope should be -2
# extra fun: prove the conjecture above
# Continue with this I guess: nevermind, there's a bit of maths there
# Tomatos genis??   

# TODO
# Do plot with wp_p (FREE)
# Look through some values of d
# For d = 1
# do the presentation
# - motivation
# - theoretical literature
# - confuse them so they think it's reeeeeally cool and ur smart