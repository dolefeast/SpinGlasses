module calculations

include("metropolis.jl")
include("IOutils.jl")

using .metropolis
using .IOutils

using AliasTables
using Plots

global nPoints = 100;
global nBonds = 100;
global sigma = 1.;
global T = 1.;

export nPoints
export nBonds;
export sigma;
export T;

function mean(arr)
	avg = sum(arr)/length(arr)
	return avg	
end

function std(arr, arrmean)
	return sqrt(mean( (arr .- arrmean).^2))
end

function arrayStatistics(arr, f=identity)
	# statistics of a f_m: N -> R, f_m(n) := arr[m][n] for the realization m at index n
	realizations = length(arr)

	maxN = length(arr[1])
	meanArr = [ 0. for _ in 1:maxN ]
	stdArr = [ 0. for _ in 1:maxN ]
	for n in 1:maxN
		mean = 0.
		std = 0.
		for realization in arr
			mean += f(realization[n])
		end
		mean /= realizations
		for realization in arr
			std += (f(realization[n]) - mean) ^2 # needs to be divided
		end
		std = sqrt(std/(realizations-1)/realizations)

		meanArr[n] = mean
		stdArr[n] = std
	end

	return meanArr, stdArr
end

function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r * sin.(θ), k .+ r * cos.(θ)
end

function addCirclePlot(h, k, r)
    """
    Plots a circle of radius r and center at (h, k)
    """
    return plot!(circleShape(h,k,r), seriestype=[:shape], lw=2,
		  c = :blue, linecolor = :black, legend=false, fillalpha = 0., aspect_ratio = 1, grid=false, ticks=false, axes=false, axis=([], false))
end # function

function pol2car(r, θ)
    r*cos(θ), r*sin(θ)
end # function

function nodesCartesian()
    """
    Draws the unit circle, and nPoints evenly distributed along it.
    Returns the list of cartesian coordinates of each of the points
    """
    points_θ = LinRange(π/2, 5π/2, nPoints + 1); # I want to start counting at the top node
    pointsCartesian = @. pol2car(1, points_θ);

    if nPoints < 20 # For visualization purposes
        scatter!(pointsCartesian)
    end
    return pointsCartesian
end

function connectionsCartesian(pointsCartesian, connectionsArray)
    """
    Returns the bonds in cartesian cordinates
    """
    lines = []
    for (startingPointIndex, endingPointIndex) in connectionsArray
        append!(lines, [[pointsCartesian[startingPointIndex], 
                         pointsCartesian[endingPointIndex]]])
    end
    lines
end

function drawCircleAndPoints(nPoints::Int64)
    """
    Draws the unit circle, and nPoints evenly distributed along it.
    Returns the list of cartesian coordinates of each of the points
    """
    p = addCirclePlot(0,0,1)
    points_θ = LinRange(π/2, 5π/2, nPoints + 1);
    pointsCartesian = @. pol2car(1, points_θ);

    if nPoints < 20
        scatter!(pointsCartesian)
    end
    return pointsCartesian
end

function errorPlot(x, y)
	mean = y[1]
	error = y[2]
	plot(x, mean, ribbon = (error))
end

function generateConnectionSet(arrayFormat::Bool=false)
    """
	Draws nBonds bonds among nPoints nodes, with probability given by r^-(1+sigma), with r the distance between the nodes. For higher dimensionality d one would have r^-(d+sigma)
    """
	
	if nBonds >div( nPoints * (nPoints - 1) , 2)
		error("nBonds > nPoints * (nPoints - 1) / 2")
	end
	
	points = 1:(nPoints - 1)
	distance = min.( points, div(nPoints, 2), nPoints .- points) 
    startingPointAT = AliasTables.AliasTable(ones(nPoints))
	distanceAliasTable = AliasTables.AliasTable((1 ./ float.(distance) .^ (1+sigma)))
    connectionsArray = Set()

    while length(connectionsArray) < nBonds
        particle1Choice = rand(startingPointAT, 1)[1]
        particle2Addition = rand(distanceAliasTable, 1)[1]
        particle2Choice = particle2Addition + particle1Choice
        if particle2Choice > nPoints
            particle2Choice = particle2Choice%nPoints 
        end
        push!(connectionsArray, Set([particle1Choice, particle2Choice ]))
    end  # while
	if arrayFormat
		return 	connectionSetToConnectivityArrayConversion(connectionsArray)
	end
    connectionsArray
end # function 

function correlationFunction(connectionsArray)
    """
    The two-point correlation function of a connection Array
    is distribution of the number of connections given for 
    a certain distance
    Parameters:
    connectionsArray = Array{Tuple{Int64,Int64}[]}[] # The syntax is messed up
    Returns
    correlationFunction = Array{Int64}[nPoints]
    """
    histogram = [0 for _ in 1:div(nPoints, 2)]
    for (startingPoint, endingPoint) in connectionsArray
		if startingPoint > endingPoint # Since connections array is stored as a set (unordered)
            c = endingPoint
            endingPoint = startingPoint
            startingPoint = c
        end
        distance = min(endingPoint - startingPoint, startingPoint + nPoints - endingPoint)  
		# println("$startingPoint and $endingPoint are separated by $distance")
        #print( "\nStarting point: ", startingPoint, ", ending point: ", endingPoint, "\n\tDistance: ", distance)
        histogram[distance] += 1
    end
    histogram ./ nBonds
end
	
function correlationFunctionStatistics(realizations)
    maxDistance = div(nPoints, 2)+1
    realizationsCorrelationFunction = [ [ 0. for _ in 1:realizations ] for _ in 1:maxDistance ]
    averageCorrelationFunction = [ 0. for _ in 1:maxDistance ]
    stdCorrelationFunction = [ 0. for _ in 1:maxDistance ]
    for realization in 1:realizations
        connect = generateConnectionSet()
        cF = correlationFunction(connect)
        for (i, value) in enumerate(cF)
            realizationsCorrelationFunction[i][realization] = value
        end
    end

    for i in 1:maxDistance
         averageCorrelationFunction[i] = mean(realizationsCorrelationFunction[i])
         stdCorrelationFunction[i] = std(realizationsCorrelationFunction[i])./sqrt(realizations)
    end
    (averageCorrelationFunction, stdCorrelationFunction)
end

function whichCluster(node, clusterArray)
    for (i, cluster) in enumerate(clusterArray)
            if node in cluster
            return i
        end
    end
    error("Node was not in the cluster") # This error should never be achieved
end

function bigPlotFunction(connectionSet, clusterColoring=false)
    """
    Plots connections randomly generated
    """
	p = plot()
	plot!(p, addCirclePlot(0, 0, 1))
    pointsCartesian = drawCircleAndPoints(nPoints)
    bondsCartesian = connectionsCartesian(pointsCartesian, connectionSet)

    if clusterColoring
        connectionArray = connectionSetToConnectivityArrayConversion(connectionSet)
        clusters = clusterIdentification(connectionArray)
        colors = [ ( "#" * join(rand("ABCDEF0123456789", 6)) ) for _ in clusters ]
    end
    for (bondPlot, bond) in zip(bondsCartesian, connectionSet)
        if clusterColoring
            node = first(bond) # since bond is a set.
            colorIndex = whichCluster(node, clusters)
            color = colors[colorIndex]
        else
            color = :purple
        end
        p = plot!(p, bondPlot, color=color, linewidth=3)
    end
	p
end # function

function clusterIdentification(connectivityArray)
    """
    The connectivityArray is an array in which connectivityArray[i] is an array that says to which particles the particle i is connected to
    clusterIdentification returns an nPoints long array for which each entry is the cluster to which the particle belongs
    The algorithm goes as follows:
        Initialization:
                clusterTag = 0
                1. Assign to every particle the cluster 0
        Procedure if the particle i has no cluster, assign to it the current clusterTag
                To each particle that is connected to the particle i assign also the cluster clusterTag


        """

    nPoints = length(connectivityArray)
    visited = Set{Int}()
    clusters = []

    function dfs(node, cluster)
        stack = [ node ]
        while !isempty(stack)
            curr = pop!(stack)
            if !(curr in visited)
                push!(visited, curr)
                push!(cluster, curr)
                for neighbor in connectivityArray[curr]
                    if !(neighbor in visited)
                        push!(stack, neighbor)
                    end
                end
            end
        end
        cluster
    end

        for i in 1:nPoints
            if !(i in visited)
                cluster = Set{Int}()
                cluster = dfs(i, cluster)
                append!(clusters, [cluster])
            end
        end
        clusters
end #clusterIdentification

function clusterSizeDistribution(clusters::Vector{Any})
    """
    Calculates the distribution of cluster sizes for a given connectivity array
    """

    clusterSizeArray = [ 0 for _ in 1:nPoints ] # The biggest possible cluster containing every point

    for cluster in clusters
        clusterSizeArray[length(cluster)] += 1
    end
    clusterSizeArray
end

function connectionSetToConnectivityArrayConversion(connectionSet)
	"""
	Converts the conectionSet into a connectivity array. The set::Set{Set{Int, Int}} contains information of which node is connected to which node. The connectivityArray is of length nPoints and each entry is the nodes to which the pointed node is connected to.
	"""
    connectivityArrayUnsorted = [ [] for _ in 1:nPoints ]
    connectivityArraySorted = [ [] for _ in 1:nPoints ]
    for (startingPoint, endingPoint) in connectionSet
        append!(connectivityArrayUnsorted[startingPoint], endingPoint)
        append!(connectivityArrayUnsorted[endingPoint], startingPoint)
       end
        
       for (i, connection) in enumerate(connectivityArrayUnsorted)
           connectivityArraySorted[i] = sort(connection)
       end

    connectivityArraySorted    
end

function clusterSizeStatistics(realizations::Int64=100)
	"""
	Calculates the mean and std of the mean of the dstributions of clusterSizes
	"""
	# The Maximum cluster size is nPoints
	lengthDistributionsRealizations = [ []  for _ in 1:realizations]

	for realization in 1:realizations
		connections = generateConnectionSet(true)	
		clusters = clusterIdentification(connections)
		sizesSet = clusterSizeDistribution(clusters)	
		
		lengthDistributionsRealizations[realization] = sizesSet
	end

	return arrayStatistics(lengthDistributionsRealizations)
end

function clusterCountStatistics(realizations::Int64=100)
	"""
	Calculates the statistics of the number of clusters
	"""
	nClusters = [ 0 for _ in 1:realizations ]	
	for realization in 1:realizations
		connections = generateConnectionSet(true)	
		clusters = clusterIdentification(connections)
		nClusters[realization] = length(clusters)
	end # for

	arrayStatistics(nClusters)
	end # function


end # module