using Statistics;

function std(arr)
    avg = mean(arr)
    sqrt( sum( (arr .- avg) .^ 2 ) / (length(arr) - 1))
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
    histogram = [0 for _ in 1:nPoints]
    for (startingPoint, endingPoint) in connectionsArray
        distance = min(abs(endingPoint-startingPoint), abs(startingPoint+nPoints-endingPoint))  
        histogram[distance] += 1
    end
    histogram
end

function averageCorrelationFunction(nPoints, maxNumberOfBonds, interactionExponent, connectionsArrayMethod, nRealizations = 1000)
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
    realizationsCorrelationFunction = [ [ 0 for _ in 1:nRealizations ] for _ in 1:nPoints ]
    averageCorrelationFunction = [ 0. for _ in 1:nPoints ]
    stdCorrelationFunction = [ 0. for _ in 1:nPoints ]
    for realization in 1:nRealizations
        connect = connectionsArrayMethod(nPoints, maxNumberOfBonds, interactionExponent)
        cF = correlationFunction(nPoints, connect)
        for (i, value) in enumerate(cF)
            realizationsCorrelationFunction[i][realization] = value
        end
    end

    for i in 1:nPoints
         averageCorrelationFunction[i] = mean(realizationsCorrelationFunction[i])
         stdCorrelationFunction[i] = std(realizationsCorrelationFunction[i])
    end
    (averageCorrelationFunction, stdCorrelationFunction)
end
