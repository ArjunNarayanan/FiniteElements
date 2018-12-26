module boundary_conditions

using assembly

export applyDirichletBCs


"""
	ApplyDirichletBCs(global_system::GlobalSystem, index::Int64,
	value::Float64)	
Modify `global_system` so that the solution degree of 
freedom corresponding to the global `index` is assigned
`value`.
"""
function applyDirichletBCs(global_system::GlobalSystem, 
					index::Int64,
					value::Float64)
	modify_rhs = global_system.K[:,index]
	for i in modify_rhs.nzind
		if (i == index)
			global_system.F[i] = modify_rhs[i]*value
		else
			global_system.F[i] -= modify_rhs[i]*value
			global_system.K[i, index] = 0.0
			global_system.K[index, i] = 0.0
		end
	end
end



"""
	ApplyDirichletBCs(system::GlobalSystem,
			boundary_node_ids::Array{Int64, 1},
			boundary_node_dofs::Array{Array{Int64, 1}, 1},
			boundary_values::Array{Array{Float64, 1}, 1})
Enforce `boundary_values` at the degrees of freedom specified by
`boundary_node_ids` and `boundary_node_dofs`.
# Example
	boundary_node_ids = [1,2]
	boundary_node_dofs = [[1],
						  [1,2]]
	boundary_values = [[0.5],
	  				   [0.1,0.8]]
enforces the value `0.5` at degree of freedom `1` of
node `1`.
"""
function applyDirichletBCs(system::GlobalSystem,
			boundary_node_ids::Array{Int64, 1},
			boundary_node_dofs::Array{Array{Int64, 1}, 1},
			boundary_values::Array{Array{Float64, 1}, 1})
	@assert length(boundary_node_ids) == length(boundary_node_dofs)
	@assert length(boundary_node_ids) == length(boundary_values)
	for i in eachindex(boundary_node_dofs)
		@assert length(boundary_node_dofs[i]) == length(boundary_values[i])
	end

	for I in eachindex(boundary_node_ids)
		nodeID = boundary_node_ids[I]
		for i in eachindex(boundary_node_dofs[I])
			dof = boundary_node_dofs[I][i]
			global_index = (nodeID - 1)*system.ndofs + dof
			value = boundary_values[I][i]
			applyDirichletBCs(system, global_index, value)
		end
	end

end



# module boundary conditions ends here
end
# module boundary conditions ends here