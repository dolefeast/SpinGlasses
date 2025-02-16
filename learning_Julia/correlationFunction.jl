using Statistics;

function std(arr)
    avg = mean(arr)
    sqrt( sum( (arr .- avg) .^ 2 ) / (length(arr) - 1))
end

function errorPlot(y, color=:black)
    yMean = y[1]
    yTop = y[1] .+ y[2] ./2
    yBot = y[1] .- y[2] ./2

    plot!(yMean, linecolor=color, legend=false)
    plot!(yTop, linecolor=color, ls=:dash, legend=false)
    plot!(yBot, linecolor=color, ls=:dash, legend=false)
end

function correlationFunction(nPoints, connectionsArray)
    """
    The two-point correlation function of a connection Array
    is distribution of the number of connections given for 
    a certain distance
    Parameters:
    connectionsArray = Array{Tuple{Int64,Int64}[]}[] # The syntax is messed up
    Returns
    correlationFunction = Array{Int64}[nPoints]
    """
    nBonds = length(connectionsArray)
    histogram = [0 for _ in 1:div(nPoints, 2)]
    for (startingPoint, endingPoint) in connectionsArray
        if startingPoint > endingPoint
            c = endingPoint
            endingPoint = startingPoint
            startingPoint = c
        end
        distance = min(endingPoint - startingPoint, startingPoint + nPoints - endingPoint)  
        #print( "\nStarting point: ", startingPoint, ", ending point: ", endingPoint, "\n\tDistance: ", distance)
        histogram[distance] += 1
    end
    histogram ./ nBonds
end

function averageCorrelationFunction(nPoints, maxNumberOfBonds, interactionExponent, connectionsArrayMethod, nRealizations = 10000)
    """
    Parameters
        nPoints: Int64, the number of points to consider in the array
        maxNumberOfBonds: Int64, the number of bonds
        interactionExponent: Float64, the exponent of the interaction 1/d**e
        connectionsArrayMethod: The method with which the connections array are realised
        nRealizations: how many times the correlation function is calculated
    Returns:
        averageSTDCorrelationFunction = [averageCorrelationFunction, stdCorrelationFunction]
    """
    maxDistance = div(nPoints, 2)
    realizationsCorrelationFunction = [ [ 0. for _ in 1:nRealizations ] for _ in 1:maxDistance ]
    averageCorrelationFunction = [ 0. for _ in 1:maxDistance ]
    stdCorrelationFunction = [ 0. for _ in 1:maxDistance ]
    for realization in 1:nRealizations
        connect = connectionsArrayMethod(nPoints, maxNumberOfBonds, interactionExponent)
        cF = correlationFunction(nPoints, connect)
        for (i, value) in enumerate(cF)
            realizationsCorrelationFunction[i][realization] = value
        end
    end

    for i in 1:maxDistance
         averageCorrelationFunction[i] = mean(realizationsCorrelationFunction[i])
         stdCorrelationFunction[i] = std(realizationsCorrelationFunction[i])
    end
    (averageCorrelationFunction, stdCorrelationFunction)
end

function connectionSetToConnectivityArrayConversion(nPoints, connectionSet)
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

