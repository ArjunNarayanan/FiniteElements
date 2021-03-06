module assembly

using geometry, SparseArrays

export Assembler, SystemMatrix, SystemRHS,
		updateSystemMatrix, updateSystemRHS,
		GlobalSystem, extract






"""
	elementMatrix(T::Triangulation{N,dim}, dofs::Int64) where {N,dim}
Return an array of size `(N,N)` (the element matrix) each of
whose entries is a `(dofs,dofs)` zero matrix.
"""
function elementMatrix(T::Type{<:Triangulation{N,dim}},
	dofs::Int64) where {N,dim}
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
function elementRHS(T::Type{<:Triangulation{N,dim}},
	dofs::Int64) where {N,dim}
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
"""
struct SystemMatrix
	I::Array{Int64, 1}
	J::Array{Int64, 1}
	vals::Array{Float64, 1}
	function SystemMatrix()
		I = Array{Int64, 1}()
		J = Array{Int64, 1}()
		vals = Array{Float64, 1}()
		new(I, J, vals)
	end
end


"""
	SystemRHS
Struct to store the information necessary to construct a sparse
vector for the right hand side.
# Fields
- `I::Array{Int64, 1}` row index
- `vals::Array{Float64, 1}` corresponding values.
"""
struct SystemRHS
	I::Array{Int64, 1}
	vals::Array{Float64, 1}
	function SystemRHS()
		I = Array{Int64, 1}()
		vals = Array{Float64, 1}()
		new(I, vals)
	end
end



"""
	updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	nodes::Array{Int64, 1}, ndofs::Int64)
Update `system_matrix` with the values in `element_matrix`.
`nodes` is an array of global node numbers. `ndofs` is the
number of degrees of freedom per node.
"""
function updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	nodes::Array{Int64, 1}, ndofs::Int64)

	for I in 1:length(nodes)
		node_I = nodes[I]
		for J in 1:length(nodes)
			node_J = nodes[J]
			counter = 1
			for j in 1:ndofs
				global_j = (node_J - 1)*ndofs + j
				for i in 1:ndofs
					global_i = (node_I - 1)*ndofs + i

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
	updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	nodes1::Array{Int64, 1}, nodes2::Array{Int64, 1})
Update `system_matrix` with the values in `element_matrix`.
`nodes1` and `nodes2` are an array of global node numbers.
The method assumes that `element_matrix[I,J]` is associated
with the coupling between global nodes `nodes1[I]` and
`nodes2[J]`. `ndofs` is the number of degrees of freedom per node.
"""
function updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	nodes1::Array{Int64, 1}, nodes2::Array{Int64, 1},
	ndofs::Int64)
	for I in 1:length(nodes1)
		node_I = nodes1[I]
		for J in 1:length(nodes2)
			node_J = nodes2[J]
			counter = 1
			for j in 1:ndofs
				global_j = (node_J - 1)*ndofs + j
				for i in 1:ndofs
					global_i = (node_I - 1)*ndofs + i

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
	nodes::Array{Int64, 1})
Update the `system_rhs` with the corresponding entries from
`element_rhs`. Use `triangulation` to get the global node
numbers. Use `system_rhs.dofs` to get the number of degrees of
freedom per node in order to compute the global dof number.
"""
function updateSystemRHS(system_rhs::SystemRHS,
	element_rhs::Array{Array{Float64}, 1},
	nodes::Array{Int64, 1}, ndofs::Int64)
	for I in 1:length(nodes)
		node_I = nodes[I]
		for i in 1:ndofs
			global_i = (node_I - 1)*ndofs + i

			value = element_rhs[I][i]

			push!(system_rhs.I, global_i)
			push!(system_rhs.vals, value)
		end
	end
end


"""
	Assembler
# Fields
	system_matrix::SystemMatrix
	system_rhs::SystemRHS
	element_matrix::Array{Array{Float64}, 2}
	element_rhs::Array{Array{Float64}, 1}
	ndofs::Int64
# Constructor
	Assembler(T::Type{<:Triangulation}, ndofs::Int64)
# Description
Collects all the structures for assembling the linear system.
"""
struct Assembler
	element_matrix::Array{Array{Float64}, 2}
	element_rhs::Array{Array{Float64}, 1}
	ndofs::Int64
	function Assembler(T::Type{<:Triangulation},
						ndofs::Int64)
		element_matrix = elementMatrix(T, ndofs)
		element_rhs = elementRHS(T, ndofs)
		new(element_matrix, element_rhs, ndofs)
	end
end



"""
	GlobalSystem
Stores the global sparse matrix `K`, the global sparse vector
of the right hand side `F`, and the global solution vector `D`.
# Attributes
	K::SparseMatrixCSC{Float64, Int64}
	D::Array{Float64, 1}
	F::SparseVector{Float64, Int64}
	ndofs::Int64
"""
mutable struct GlobalSystem
	K::SparseMatrixCSC{Float64, Int64}
	D::Array{Float64, 1}
	F::SparseVector{Float64, Int64}
	ndofs::Int64
	"""
		GlobalSystem(system_matrix::SystemMatrix, system_rhs::SystemRHS,
			dofs_per_node::Int64, number_of_nodes::Int64)
	assemble the global linear system from the `system_matrix`, `system_rhs`.
	"""
	function GlobalSystem(system_matrix::SystemMatrix, system_rhs::SystemRHS,
		dofs_per_node::Int64, number_of_nodes::Int64)

		ndofs = dofs_per_node*number_of_nodes
		K = sparse(system_matrix.I, system_matrix.J, system_matrix.vals, ndofs,
			ndofs)
		F = sparsevec(system_rhs.I, system_rhs.vals, ndofs)
		D = zeros(F.n)
		new(K, D, F, dofs_per_node)
	end
end


"""
	extract(system::GlobalSystem,
			dof::Int64,
			node_ids::Array{Int64, 1})
Extract the solution values from `system.D` corresponding to the `dof` degree of freedom,
and the node numbers corresponding to `node_ids`.
"""
function extract(system::GlobalSystem, dof::Int64,
				 node_ids::Array{Int64, 1})
	indices = [((I-1)*system.ndofs + dof) for I in node_ids]
	return system.D[indices]
end

"""
	extract(system::GlobalSystem, dof::Int64,
				 node_id::Int64)
Extract the solution value from `system.D` corresponding
to the `dof` degree-of-freedom and `node_id` node.
"""
function extract(system::GlobalSystem, dof::Int64,
				 node_id::Int64)
	index = (node_id - 1)*system.ndofs + dof
	return system.D[index]
end


"""
	extract(system::GlobalSystem, node_id::Int64)
Extract the solution value from `system.D` corresponding
to all degrees of freedom of node `node_id`.
"""
function extract(system::GlobalSystem, node_id::Int64)
	indices = [((node_id-1)*system.ndofs + I) for I in 1:system.ndofs]
	return system.D[indices]
end






# module assembly ends here
end
# module assembly ends here
