using StatsBase

function costMatrixjl(muhat, nuhat, p, d)
    mhsize = size(muhat)[1]
    nhsize = size(nuhat)[1]
    
    println("mhsize: ", mhsize)
    println("nhsize: ", nhsize)
    out = zeros(mhsize, nhsize)
    for i in 1:mhsize
        for j in 1:nhsize
            out[i,j] = (mean(abs.(muhat[:, i] - nuhat[:, j]).^q))^(p / q)
        end
    end
    return out
end