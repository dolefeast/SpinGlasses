using Plots;

include("dictionaryPlot.jl");
include("correlationFunction.jl");
include("clusterCounting.jl");

function clusteringStatistics(nPoints, nBonds, interactionExponent, connectivityMethod, nRealizations = 2000)
    nClusters = [ 0 for _ in 1:nRealizations ]
    for i in 1:nRealizations
        # A set of the different bonds. The length of this set is nBonds.
        connectiononivityMethod(nPoints, nBonds, interactionExponent)
        # An array of length nPoints. At its entry n, another array (min length 0, max length nPoints(nPoints - 1)/2), with every neighbor of the point n. 
        connectionArray = connectionSetToConnectivityArrayConversion(nPoints, connectionSet)
        nClusters[i] = clusterIdentification(connectionArray)     
    end 
    mean(nClusters), std(nClusters)
end
