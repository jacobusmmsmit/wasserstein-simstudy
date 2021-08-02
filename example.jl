using OptimalTransport
using Distances
using Tulip
using Distributions

# have 200 pts, give them equal weight

M = 200
μ = fill(1 / M, M)
μsupport = rand(M)

N = 250
ν = fill(1 / N, N)
νsupport = rand(N);

C = pairwise(SqEuclidean(), μsupport', νsupport'; dims=2);

μ
ν
νsupport

N = 10
μ = Normal(0, 1)
ν = Normal(0, 1)
μ̂ = sort(rand(μ, N))
ν̂ = sort(rand(ν, N))

# Estimated but quick optimal (shortcut known):
√(mean((μ̂ - ν̂).^2))

# Estimated using no shortcuts:
v2m(sample) = reshape(sort(sample), 1, length(sample))
C = pairwise(SqEuclidean(), v2m(μ̂), v2m(ν̂); dims=2)

√(emd2(
    fill(1 / N, N),
    fill(1 / N, N), 
    C, 
    Tulip.Optimizer()
))

using BenchmarkTools
@benchmark √(emd2(fill(1 / N, N), fill(1 / N, N), C, Tulip.Optimizer()))

# Investigate:
# Easy ("Good"): Look at variance of Wp^p (wasserstein p distance to the p) and Wp when μ = ν as N increases
# Loop, store distance, compute estimator of sample variance (1/n × sum((xi - x̄)^2))


# , there are CLT for cost 



# using AssignmentSolver

# cost = C
# reward = Float64.(reward2cost(cost))
# sqrt(compute_objective(hungarian_assignment(reward), cost))