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
        length(clusters)
end
