set.seed(1)
library(tidyverse)
library(glue)
library(Rcpp)
library(RcppEigen)
library(scales)
sourceCpp("src/network_simplex_fast.cpp") # Compile and load the solver. The header files are also necessary.
sourceCpp("src/squared_cost_eval.cpp")
# sourceCpp("src/costMatrix.cpp")

# n <- 100  # Sample size (unused)
d <- 1    # Dimension
p <- 2L    # Wasserstein-p norm
q <- 2L    # L-q metric to calculate cost matrix
K <- 2000 # Number of replicates
R <- 1000 # Number of bootstrap replicates
ns <- round(10^(seq(2, 3.5, 1/4)))
# ns <- c(100, 300, 1000, 3000) #, 10000)
ds <- c(1, 10, 100)   # Dimension
var_wp_emp <- rep(NA, length(ns))
var_wp_emp_p <- rep(NA, length(ns))
wp_boot_sd <- rep(NA, length(ns))
wp_boot_p_sd <- rep(NA, length(ns))
results <- tibble()
for (d in c(7)){
    print(glue("Dimension: ", d))
    j <- 0
    for (n in ns){
        print(glue("No. samples: ", n))
        wp_emp_p <- rep(NA, K)
        j <- j+1
        for (k in 1:K){
            # 1.) Generate samples
            muhat <- matrix(rnorm(n * d), nrow = d)
            nuhat <- matrix(rnorm(n * d), nrow = d)
            
            # wp_emp_p sort, diff, mean, it's somewhere bro just like find it
            
            # 2.) Compute the cost matrix 
            cost_mat <- EvaluateSquaredCost(muhat, nuhat, 4L) # d(x, y) = L-q metric, Wasserstein-p norm # costMatrix(muhat, nuhat, q, p)
            # cost_mat <- costMatrix(muhat, nuhat, q, p)
            
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
    
    results <- bind_rows(results, tibble(d = d, ns, var_wp_emp, wp_emp_lci, wp_emp_uci))
    coefs <- lm(log10(var_wp_emp) ~ log10(ns))$coefficients
    
    # qplot(ns, var_wp_emp) +
    #     theme_bw() +
    #     scale_x_log10() +
    #     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #                   labels = trans_format("log10", math_format(10^.x))) +
    #     geom_ribbon(aes(x = ns, ymin = wp_emp_lci, ymax = wp_emp_uci), alpha = 0.3, fill = "#008ccc") +
    #     geom_abline(intercept = coefs[1], slope = coefs[2], linetype = "dashed", size = 0.7) +
    #     labs(x = "Number of Samples",
    #          y = "Variance",
    #          title=glue("Variance of Wasserstein-2 Estimator in ", d," Dimensions"),
    #          subtitle = glue("Absolute value of gradient = ", as.character(round(-coefs[2], 3))))
    # 
    # output_dir = "outputs/"
    # save_plot_as = glue("wasserstein2_", d, "_plot.pdf")
    # ggsave(glue(output_dir, save_plot_as), device = "pdf", width = 5, height = 4)
}
models <- results %>%
    group_by(d) %>%
    do(model = lm(log10(var_wp_emp) ~ log10(ns), data = .)$coefficients[[2]]) %>%
    mutate(c = unlist(model)) %>%
    select(-model)

models %>%
    ggplot(aes(d, -c)) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 1.2) +
        scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
        theme_bw() +
        labs(x = "Dimension of underlying spaces",
             y = substitute(a*b*c, list(a = "c (Variance decays at rate  ", b = quote(frac(1, n^c)), c = " )")),
             title = "Rate of convergence of variance")
ggsave(glue(output_dir, "final_plot.pdf"), device = "pdf", width = 5, height = 4)