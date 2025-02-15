function clusterIdentification(connectivityArray)
    """
    The connectivityArray is an array in which connectivityAray[i] is an array that says to which particles the particle i is connected to
    clusterIdentification returns an nPoints long array for which each entry is the cluster to which the particle belongs
    The algorithm goes as follows:
        Initialization:
                clusterTag = 0
                1. Assign to every particle the cluster 0
        Procedure
                if the particle i has no cluster, assign to it the current clusterTag
                To each particle that is connected to the particle i assign also the cluster clusterTag


        """

    nPoints = length(connectivityArray)
    clusterTag = 1
    clustering = [ 0 for _ in 1:nPoints ]
       
    change = false
    for (i,outgoingArray) in enumerate(connectivityArray)
        if clustering[i] == 0
            clustering[i] = clusterTag
            clusterTag += 1
        end
        for (j, connection) in enumerate(outgoingArray)
            if clustering[connection] != 0
            end
            clustering[connection] = clustering[i]

        end
        print(i, ", ", clustering, "\n")
    end
    clustering
end

