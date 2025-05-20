module  metropolis

function singleStep(nPoints::Int64, T::Float64, connections::Vector{Vector{Any}}, spins::Vector{Int64})
	"""
	Computes single Metropolis step based on the "long range" Ising model. 
	"""

	p = rand(1:nPoints) # Selects the particle whose spin might be flipped
	spin = spins[p] # the  spin of  the particle p
	totalE = 0.
	for neighbor in connections[p]
		totalE += spins[neighbor]
	end

	totalE *= spin

	if totalE < 0
		spins[p] *= -1
	else
		if rand() < exp(-totalE / T)
			spins[p] *= -1
		end
	end
	spins
end

function evolve(connections::Vector{Vector{Int64}}, spins::Vector{Int64}, realizations=500)
	"""
	Calculates the Metropolis time evolution of a system. It is basically Ising model with a different 
	topology.
	Parameters
		connections::Vector{Vector{Int64}}, the connections array (site [i] is the connections of the node i)
		spins::Vector{Int64}, the starting spins array
	"""
	if length(spins) != nPoints
		error("length(spins) != nPoints")
	end # if

end # function
end # module

