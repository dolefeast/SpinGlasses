using AliasTables
"""
How to plot a fixed number of connections: 
        Using Walker's alias method, draw a random point. From this random point, choose another random point. If that connection has already been done, dismiss it. 
    """
function distancesArray(nPoints)
    """
    Parameters:
        nPoints: int. The number of points in the system
    Returns:
        distances: An Array{Float}() object encoding the distances to each other point in the system. Kind of encases the topology.
    """
    distances = []
    nPoints = nPoints
    distances = 1:div(nPoints,2)
    distances
end

function pregenerateBondsList(nPoints, maxNumberOfBonds)
    bondsListAT = AliasTable(ones(nPoints))
    sampling = rand(bondsListAT, maxNumberOfBonds)
    bondsList = [ 0 for _ in 1:nPoints ]
    for bond in sampling
        bondsList[bond] += 1
    end
    bondsList
end

function walkerAliasPregeneratedConnectionsArray(nPoints, maxNumberOfBonds, interactionExponent)
    # the aliasTable gives the probability of a certain point of bonding with a neighbout.
    # it should therefore be of nPoints-1 entries.
    # My algorithm randomly chooses a point, and then assigns to this point an addition index, so that the addition is given by the interaction exponent.
    #
    # 1. Generate a bondsList: How many bonds each node is going to get
    # 2. Go through each node and generate as many bonds as the bondsList says
    #           2.1. These bonds only go "forward"i.e. to the 
    #           next nPoints/2 nodes
    # 3. Before each bonding, reinitialize the AliasTable, so that the probability of bonding with an already bonded node is 0.
    #
    # An nPoints-1 long array
    
    # An nPoints long array which contains the information of how many
    # bonds each node has
    nOfBondsList = pregenerateBondsList(nPoints, maxNumberOfBonds)
    bondsList = []

    distanceAliasTable = AliasTable(1 ./ distancesArray(nPoints) .^ interactionExponent)
       
    for (particle1, nBonds) in enumerate(nOfBondsList)
        if nBonds == 0
            continue
        end
        particle1BondsList = []
        while length(particle1BondsList) < nBonds
            particle2 = rand(distanceAliasTable, 1)[1] + particle1
            if particle2 > nPoints
                particle2 = particle2%nPoints
            end
            if particle2 in particle1BondsList
                continue
            end
            append!(bondsList, [ sort([ particle1, particle2 ]) ])
            append!(particle1BondsList, (particle2))
        end
    end

    bondsList
end

function walkerAliasConnectionsSet(nPoints, maxNumberOfBonds, interactionExponent)
    # the aliasTable gives the probability of a certain point of bonding with a neighbout.
    # it should therefore be of nPoints-1 entries.
    # My algorithm randomly chooses a point, and then assigns to this point an addition index, so that the addition is given by the interaction exponent
    #
    # An nPoints-1 long array
    chooseAT = AliasTable(ones(nPoints))
    distanceAliasTable = AliasTable(1 ./ distancesArray(nPoints) .^ interactionExponent)
    connectionsArray = Set()
    while length(connectionsArray) < maxNumberOfBonds
        particle1Choice = rand(chooseAT, 1)[1]
        particle2Addition = rand(distanceAliasTable, 1)[1]
        particle2Choice = particle2Addition + particle1Choice
        if particle2Choice > nPoints
            particle2Choice = particle2Choice%nPoints 
        end
        # Use x*nPoints + y as a key of the lookup dictionary, where 
        #       x = max(particle1Choice, particle2Choice)
        #       y = min(particle1Choice, particle2Choice)
        push!(connectionsArray, Set([particle1Choice, particle2Choice ]))
    end
    connectionsArray
end


function walkerAliasConnectionsArray(nPoints, maxNumberOfBonds, interactionExponent)
    # the aliasTable gives the probability of a certain point of bonding with a neighbout.
    # it should therefore be of nPoints-1 entries.
    # My algorithm randomly chooses a point, and then assigns to this point an addition index, so that the addition is given by the interaction exponent
    #
    # An nPoints-1 long array
    chooseAT = AliasTable(ones(nPoints))
    distanceAliasTable = AliasTable(1 ./ distancesArray(nPoints) .^ interactionExponent)
    connectionsArray = []
    while length(connectionsArray) < maxNumberOfBonds
        particle1Choice = rand(chooseAT, 1)[1]
        particle2Addition = rand(distanceAliasTable, 1)[1]
        particle2Choice = particle2Addition + particle1Choice
        if particle2Choice > nPoints
            particle2Choice = particle2Choice%nPoints 
        end
        # Use x*nPoints + y as a key of the lookup dictionary, where 
        #       x = max(particle1Choice, particle2Choice)
        #       y = min(particle1Choice, particle2Choice)
        if (particle1Choice, particle2Choice) in connectionsArray
            continue
        end
        append!(connectionsArray, [(particle1Choice, particle2Choice) ])
    end
    connectionsArray
end

function firstElement(arr)
    arr[1]
end

function identifyGroups(nPoints,connectionsArray)
    """
    Given a connections array identify each particle with a group.
    Parameters:
        nPoints: the amount of particles in the system
        connectionsArray: Which particles are connected to which
    Returns:
        group # an array of integers of length nPoints so that group[i] is the group to which the particle i belongs
    
    Both nPoints and connectionsArray must be provided as connectionsArray does not encode how many particles there are in the system. e.g. consider the system constituted by 5000 particles and just 1 bond.
    The algorithm goes: 
        1. Assign each particle a 0. Zero means no belonging to any group
        2. The first particle gets assigned a 1. Every particle connected to this particle is also assigned a 1.
        3. Go to particle 2. If this particle has no asigned group, asign it group 2, otherwise, leave as is. Asign to every particle connected to particle 2 its group.
    """
    connectionsArray = sort(connectionsArray, by=firstElement)

    groups  = zeros(Int32, nPoints)
    groupCount = 1
    unassignedParticles = 1:nPoints
    for (startingPoint, endingPoint) in connectionsArray
        if groups[startingPoint] == 0 && groups[endingPoint] == 0
            groups[startingPoint] = groupCount
            groups[endingPoint] = groupCount
            groupCount += 1
        elseif groups[startingPoint] != 0 && groups[endingPoint] == 0
            groups[endingPoint] = groups[startingPoint]
        elseif groups[startingPoint] == 0 && groups[endingPoint] != 0
            groups[startingPoint] = groups[endingPoint] 
        elseif groups[startingPoint] != 0 && groups[endingPoint] != 0
            print("Starting point ", startingPoint, ", group: ", groups[startingPoint], "\n")
            print("ending point ", endingPoint, ", group: ", groups[endingPoint], "\n")
            if groups[startingPoint] !=  groups[endingPoint]
                print("The groupings so far: ", groups)
                throw(DomainError)
            end
        end
    end
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

