# Pkg.activate("/home/sanz/spinGlasses/calculations/")

using .calculations
using CSV
using DataFrames

function  saveClusterSize(params, nRealizations)
    calculations.T, calculations.nPoints, calculations.nBonds, calculations.sigma = params

    clusterSizes = []
    for _ in 1:nRealizations
        J = calculations.generateConnectionSet()
        connectionsArray = calculations.connectionSetToConnectivityArrayConversion(J)

        clusters = calculations.clusterIdentification(connectionsArray)
        sizesSet = calculations.clusterSizeDistribution(clusters)
        push!(clusterSizes, sizesSet)
    end

    CSV.write("/home/dein/spinGlasses/juliaScripts/results/Cluster Counting/$(calculations.filename(params...)).csv", DataFrame(clusterSizes, :auto), header=false)
end


TArray = [1,]
nPointsArray = [1024, 4096, 8192, 16384]
coordNumberArray = [0.5, 0.75, 1, 1.5, 2, 4, 8, 10]
sigmaArray = [0.2, 0.8, 1.5]

nRealizations = 100

paramsArray = calculations.buildParamsArray(TArray, nPointsArray, coordNumberArray, sigmaArray)

# Average number of clusters
for params in paramsArray
    saveClusterSize(params, nRealizations)
    println(calculations.filename(params...))
end

