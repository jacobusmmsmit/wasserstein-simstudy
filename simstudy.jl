using OptimalTransport
using Distances
using Tulip
using Distributions


# uniform histograms
μ = fill(1 / 5, 5)
ν = fill(1 / 200, 200)

# random cost matrix
C = pairwise(SqEuclidean(), rand(1, 5), rand(1, 200); dims=2)

# regularization parameter
ε = 0.01

# solve entropically regularized optimal transport problem
emd2(μ, ν, C, Tulip.Optimizer())

wasserstein(Normal(2, 1), Normal(4, 1))


v2m(sample) = reshape(sort(sample), 1, length(sample))
m1 = v2m(samples_1)
m2 = v2m(samples_2)



Cn = pairwise(SqEuclidean(), m1, m2; dims=2)
emd2(samples_1, samples_2, Cn, Tulip.Optimizer())

function empirical_wasserstein(sample1, sample2; p = 2)
    sample_1 = v2m(sample1)
    sample_2 = v2m(sample2)
    println(sample_1)
    C = pairwise(SqEuclidean(), sample_1, sample_2; dims=2)
    emd2(sample_1, sample_1, C, Tulip.Optimizer())
    # return (sum((sample1 .- sample2).^p))^(1/p), emd
end

empirical_wasserstein([0.5, 0.2, 0.3], [0.2, 0.4, 0.4])

N = 10
distr_1 = Normal(0, 1)
distr_2 = Exponential(3)
samples_1 = rand(distr_1, N)
samples_2 = rand(distr_2, N)
empirical_dist = zeros(N)

for i in 1:N
    empirical_dist[i] = empirical_wasserstein(samples_1[1:i], samples_2[1:i])
end
true_dist = wasserstein(distr_1, distr_2)
p = plot(1:N, empirical_dist, label = "Empirical Distance")
hline!(p, [true_dist], label = "True Distance", ls = :dash, lc = :red)

