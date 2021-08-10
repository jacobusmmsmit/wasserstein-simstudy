/*
This file contains functions used to compute the squared Euclidean cost matrix between two empirical distributions.
    - Empirical distributions are input as matrices.
    - Within each empirical distribution, each sample is either a column vector of size #{dimensions}.
    - Input should be in column-major order #{dimensions}-by-#{samples}.
    - There is support for parallel evalaution using OpenMP. Useful for large problems.
*/

#include <RcppEigen.h>
#include <omp.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11, openmp)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;

// Version that takes the matrices X and Y with samples stored in column-major order (so d-by-n matrices)
//[[Rcpp::export(rng = false)]]
ArrayXXd EvaluateSquaredCost(const Map<MatrixXd> &x, const Map<MatrixXd> &y, const int &nthreads = 1) 
{
    ArrayXXd cost(x.cols(), y.cols());

    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for
        for (int i = 0; i < y.cols(); i++)
        {
            cost.col(i) = (x.colwise() - y.col(i)).colwise().squaredNorm();
        }
    }
    
    return cost;
}
