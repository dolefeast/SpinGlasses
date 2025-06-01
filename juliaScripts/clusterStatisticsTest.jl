using Pkg

Pkg.activate("calculations")

using calculations
using Plots

calculations.nPoints = 10
calculations.nBonds = 10
calculations.sigma = 1.0


distributions = []
for _ in 1:200
	_g = calculations.generateConnectionSet(true)
	_clusters = calculations.clusterIdentification(_g)
	_clusterSizes = calculations.clusterSizeDistribution(_clusters)
	append!(distributions, [_clusterSizes])
end

stats = calculations.arrayStatistics(distributions)

p1 = plot(distributions)
p2 = plot(stats[1], ribbon=stats[2]/2)

display((p1, p2))
