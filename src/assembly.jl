module assembly

using geometry, StaticArrays

export SystemMatrix, SystemRHS, elementMatrix,
		elementRHS, updateSystemMatrix,
		updateSystemRHS






"""
	elementMatrix(T::Triangulation{N,dim,spacedim}, dofs::Int64) where {N,dim,spacedim}
Return an array of size `(N,N)` (the element matrix) each of
whose entries is a `(dofs,dofs)` zero matrix.
"""
function elementMatrix(T::Type{<:Triangulation{N,dim,spacedim}}, 
	dofs::Int64) where {N,dim,spacedim}
	element_matrix = Array{Array{Float64}, 2}(undef, N, N)
	for i in 1:N
		for j in 1:N
			element_matrix[i,j] = zeros(dofs, dofs)
		end
	end
	return element_matrix
end

"""
	elementRHS(T::Type{<:Triangulation{N,dim,spacedim}},
	rank::Int64) where {N,dim,spacedim}
Return an array of size `(N,)` each of whose entries is 
a zero vector of length `dofs`.
"""
function elementRHS(T::Type{<:Triangulation{N,dim,spacedim}}, 
	dofs::Int64) where {N,dim,spacedim}
	element_rhs = Array{Array{Float64}, 1}(undef, N)
	for i in 1:N
		element_rhs[i] = zeros(dofs)
	end
	return element_rhs
end






"""
	SystemMatrix
Struct to store the information necessary to construct a sparse 
matrix for the given system of equations.
# Fields
- `I::Array{Int64, 1}` row index
- `J::Array{Int64, 1}` column index
- `vals::Array{Float64, 1}` corresponding values
- `dofs::Int64` number of degrees of freedom (dof) per node. 
Used to obtain the global dof number from the local dof number.
"""
struct SystemMatrix
	I::Array{Int64, 1}
	J::Array{Int64, 1}
	vals::Array{Float64, 1}
	dofs::Int64
	function SystemMatrix(dofs::Int64)
		I = Array{Int64, 1}()
		J = Array{Int64, 1}()
		vals = Array{Float64, 1}()
		new(I, J, vals, dofs)
	end
end


"""
	SystemRHS
Struct to store the information necessary to construct a sparse
vector for the right hand side.
# Fields
- `I::Array{Int64, 1}` row index
- `vals::Array{Float64, 1}` corresponding values
- `dofs::Int64` number of degrees of freedom (dof) per node. 
Used to obtain the global dof number from the local dof number.
"""
struct SystemRHS
	I::Array{Int64, 1}
	vals::Array{Float64, 1}
	dofs::Int64
	function SystemRHS(dofs::Int64)
		I = Array{Int64, 1}()
		vals = Array{Float64, 1}()
		new(I, vals, dofs)
	end
end



"""
	updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2}, 
	triangulation::Triangulation)
Update the `system_matrix` with the corresponding entries from 
the `element_matrix`. Use `triangulation` to get the global
node numbers. Use `system_matrix.dofs` to get the number of
degrees of freedom per node in order to compute the global
dof number.
"""
function updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	triangulation::Type{<:Triangulation{N,dim,spacedim}}) where {N,dim,spacedim}
	for I in 1:N
		node_I = triangulation.nodes[I]
		for J in 1:N
			node_J = triangulation.nodes[J]
			counter = 1
			for i in 1:system_matrix.dofs
				global_i = (node_I - 1)*system_matrix.dofs + i
				for j in 1:system_matrix.dofs
					global_j = (node_J - 1)*system_matrix.dofs + j
					value = element_matrix[I,J][counter]
					counter += 1

					push!(system_matrix.I, global_i)
					push!(system_matrix.J, global_j)
					push!(system_matrix.vals, value)
				end
			end	
		end
	end
end

"""
	updateSystemRHS(system_rhs::SystemRHS,
	element_rhs::Array{Array{Float64}, 1},
	triangulation::Type{<:Triangulation{N,dim,spacedim}}) where {N,dim,spacedim}
Update the `system_rhs` with the corresponding entries from 
`element_rhs`. Use `triangulation` to get the global node 
numbers. Use `system_rhs.dofs` to get the number of degrees of 
freedom per node in order to compute the global dof number.
"""
function updateSystemRHS(system_rhs::SystemRHS,
	element_rhs::Array{Array{Float64}, 1},
	triangulation::Type{<:Triangulation{N,dim,spacedim}}) where {N,dim,spacedim}
	for I in 1:N
		node_I = triangulation.nodes[I]
		for i in 1:system_rhs.dofs
			global_i = (node_I - 1)*system_rhs.dofs + i

			value = element_rhs[I][i]

			push!(system_rhs.I, global_i)
			push!(system_rhs.vals, value)
		end
	end
end




# module assembly ends here
end
# module assembly ends here