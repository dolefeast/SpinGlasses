using AliasTables
"""
How to plot a fixed number of connections: 
        Using Walker's alias method, draw a random point. From this random point, choose another random point. If that connection has already been done, dismiss it. 
"""
function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r * sin.(θ), k .+ r * cos.(θ)
end

function plotCircle(h, k, r)
    plot!(circleShape(h,k,r), seriestype=[:shape], lw=0.5,
         c = :blue, linecolor = :black, legend=false, fillalpha = 0., aspect_ratio = 1, grid=false, ticks=false, axes=false)
end

function pol2car(r, θ)
    r*cos(θ), r*sin(θ)
end

function drawCircleAndPoints(nPoints)
    """
    Draws the unit circle, and nPoints evenly distributed along it.
    Returns the list of cartesian coordinates of each of the points
    """
    plotCircle(0,0,1)
    points_θ = LinRange(π/2, 5π/2, nPoints + 1);
    pointsCartesian = @. pol2car(1, points_θ);

    if nPoints < 20
        scatter!(pointsCartesian)
    end
    return pointsCartesian
end

function distancesArray(nPoints)
    """
    Parameters:
        nPoints: int. The number of points in the system
    Returns:
        distances: An Array{Float}() object encoding the distances to each other point in the system. Kind of encases the topology.
    """
    distances = []
    nPoints = nPoints-1
    distances = 1:div(nPoints,2)
end

function pregenerateBondsList(nPoints, maxNumberOfBonds)
    bondsListAT = AliasTable(ones(nPoints))
    sampling = rand(bondsListAT, maxNumberOfBonds)
    bondsList = zeros(nPoints)
    for bond in sampling
        bondsList[bond] += 1
    end
    bondsList
end

function walkerAliasConnectionsArray(nPoints, maxNumberOfBonds, interactionExponent)
    # the aliasTable gives the probability of a certain point of bonding with a neighbout.
    # it should therefore be of nPoints-1 entries.
    # My algorithm randomly chooses a point, and then assigns to this point an addition index, so that the addition is given by the interaction exponent
    #
    # An nPoints-1 long array
    chooseAT = AliasTable(ones(nPoints))
    distanceAliasTable = AliasTable(1 ./ distancesArray(nPoints) .^ interactionExponent)
    bonds = []
    while length(bonds) < maxNumberOfBonds
        particle1Choice = rand(chooseAT, 1)[1]
        particle2Addition = rand(distanceAliasTable, 1)[1]
        particle2Choice = particle2Addition + particle1Choice
        if particle2Choice > nPoints
            particle2Choice = particle2Choice%nPoints 
        end
        if (particle1Choice, particle2Choice) in bonds
            continue
        end
        append!(bonds, [(particle1Choice, particle2Choice) ])
    end
    bonds
end

function probabilityConnectionsArray(nPoints, p=1)
    """
    outputs two arrays of the same length, with the length being the number of lines to be plotted. Right now it just joins each two points with a probability 1. This probability should be changed. 
    """
    connectingLines = []

    for i in 1:nPoints
        for j in 1:i-1
            if rand() < p
                append!(connectingLines, [(i, j)])
            end
         end
     end
     connectingLines
end

function plotConnections(pointsCartesian, connectionsArray)
    """
    pointsCartesian is a 2 x nPoints list with each element being the cartesian coordinates of each of the points.
    """
    lines = []
    for (startingPointIndex, endingPointIndex) in connectionsArray
        append!(lines, [[pointsCartesian[startingPointIndex], 
                         pointsCartesian[endingPointIndex]]])
    end
    lines
end
