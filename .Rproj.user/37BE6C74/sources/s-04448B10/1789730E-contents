# Investigate:
# Easy ("Good"): Look at variance of Wp^p (wasserstein p distance to the p) and Wp when ?? = ?? as N increases
# Loop, store distance, compute estimator of sample variance (1/n ? sum((xi - x??)^2))

# var(c(10 9 8 7 6 5 4 3 2 1))
# pattern: want plot variance against N, definintely for Wp^p goes down faster than 1/n, interesting: find out the exponent.
# it's n^-c c> 1 but how big is this and how does it depend on dimension
# How to see what size of c is:
# Imagine I know a-priori it's logvar = logc * c log(n)
# log log plot of variance against sample size then slope = c < -1
# Then take p^th root then calculate variance and then loglog plot of variance against sample size
# Pattern should be the same for whatever mertic and exponent p > 1
# Which exact metic to take and which exponent tot ake is pretty wide
# Start with euclidean norm exponent 2, then integer exponents
# Sample size up to 10,000
# Replicates: A couple

# Wasserstein takes a metric d and then (integral( d(x - y) ^p)^1/p

# d <- function(x, y, q = 2) {(sum(x - y)^q)^(1/q)} # = Lq-norm which then induces a metric but yeah

# preliminary investigation: CLT papers on Optimal transport ?? if I would like (only, otherwise do not, i repeat do not unless you would like)
# Came out february, good papers -> you must enjoy them. contains refernece to previous work which is pretty good
# https://arxiv.org/pdf/2102.06379.pdf

require(Rcpp)
sourceCpp("src/network_simplex_fast.cpp") # Compile and load the solver. The header files are also necessary.
sourceCpp("src/costMatrix.cpp")

n <- 1000 # Sample size
d <- 10   # Dimension

# 1.) Generate samples (they don't have to be normal)

muhat <- matrix(rnorm(n * d), nrow = d) # n samples from a standard normal, stored as one sample per column
nuhat <- matrix(rnorm(n * d), nrow = d)

# 2.) Compute the cost matrix 
# If x_i is the i-th sample from mu, its entries are c_{ij} = 1/n * ||x_i - y_i||_2^2.
# In general, you'll do c_{ij} = 1/n * d(x_i, y_i)^p.

cost_mat <- costMatrix(muhat, nuhat, 2, 2) # L-2 metric, wasserstein 2 norm

# remake the cost matrix

# 3.) Compute the distance from the cost matrix
w2_empirical <- sqrt(SolveAssignmentNetworkflow(cost_mat))

# Testing

  
