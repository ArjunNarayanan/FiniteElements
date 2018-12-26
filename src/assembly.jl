module assembly

using geometry, SparseArrays

export Assembler, updateSystem, GlobalSystem






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
	nodes::Array{Int64, 1})
Update the `system_matrix` with the corresponding entries from 
the `element_matrix`. Use `triangulation` to get the global
node numbers. Use `system_matrix.dofs` to get the number of
degrees of freedom per node in order to compute the global
dof number.
"""
function updateSystemMatrix(system_matrix::SystemMatrix,
	element_matrix::Array{Array{Float64}, 2},
	nodes::Array{Int64, 1}, ndofs::Int64)
	for I in 1:length(nodes)
		node_I = nodes[I]
		for J in 1:length(nodes)
			node_J = nodes[J]
			counter = 1
			for i in 1:ndofs
				global_i = (node_I - 1)*ndofs + i
				for j in 1:ndofs
					global_j = (node_J - 1)*ndofs + j
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
	system_matrix::SystemMatrix
	system_rhs::SystemRHS
	element_matrix::Array{Array{Float64}, 2}
	element_rhs::Array{Array{Float64}, 1}
	ndofs::Int64
	function Assembler(T::Type{<:Triangulation}, ndofs::Int64)
		system_matrix = SystemMatrix()
		system_rhs = SystemRHS()
		element_matrix = elementMatrix(T, ndofs)
		element_rhs = elementRHS(T, ndofs)
		new(system_matrix, system_rhs, element_matrix, element_rhs, ndofs)
	end
end




"""
	updateSystem(assembler::Assembler, nodes::Array{Int64, 1})
Calls `updateSystemMatrix` and `updateSystemRHS`.
"""
function updateSystem(assembler::Assembler, 
	nodes::Array{Int64, 1})
	updateSystemMatrix(assembler.system_matrix, 
		assembler.element_matrix,
		nodes, assembler.ndofs)
	updateSystemRHS(assembler.system_rhs, assembler.element_rhs,
		nodes, assembler.ndofs)
end



"""
	GlobalSystem
Stores the global sparse matrix `K`, the global sparse vector 
of the right hand side `F`, and the global solution vector `D`.
# Attributes
	K::SparseMatrixCSC{Float64, Int64}
	D::Array{Float64, 1}
	F::SparseVector{Float64, Int64}
"""
struct GlobalSystem
	K::SparseMatrixCSC{Float64, Int64}
	D::Array{Float64, 1}
	F::SparseVector{Float64, Int64}
	function GlobalSystem(system_matrix::SystemMatrix, 
		system_rhs::SystemRHS)
		K = sparse(system_matrix.I, system_matrix.J, system_matrix.vals)
		F = sparsevec(system_rhs.I, system_rhs.vals)
		D = zeros(F.n)
		new(K, D, F)
	end
end	




# module assembly ends here
end
# module assembly ends here