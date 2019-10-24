module boundary_conditions

using assembly

export applyDirichletBCs, Constraint, applyConstraints


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



"""
	Constraint
Represents a linear constraint equation that can be applied to a global
linear system of equations. The constraint equation is to be of
the form:
	coeff*u[index] = A1*u[i1] + A2*u[i2] + ... + inhomogeneity
`index` cannot appear on the right-hand-side.
# Attributes:
- `index::Int64` - the global index of the degree-of-freedom being
constrained.
- `coeff::Float64` - the coefficient of the degree-of-freedom being
constrained.
- `constrain_to::Array{Tuple{Int64, Float64}}` - Array of tuples
of the indices of the degrees of freedom appearing on the right-hand-side
of the constraint equation along with their coefficients.
- `inhomogeneity` - the inhomogeneous term in the constraint equation.
"""
struct Constraint
	index::Int64
	coeff::Float64
	constrain_to::Array{Tuple{Int64, Float64}}
	inhomogeneity::Float64
	function Constraint(index::Int64,
						coeff::Float64,
						constrain_to::Array{Tuple{Int64, Float64}},
						inhomogeneity::Float64)
		for term in constrain_to
			@assert index != term[1] "Index of DOF being constrained appears on right-hand-side."
		end
		new(index, coeff, constrain_to, inhomogeneity)
	end
end



"""
	applyConstraints(system::GlobalSystem,
						constraint::Constraint)
Modify `system` to ensure that the constraint equation defined by `constraint` is met.
"""
function applyConstraints(system::GlobalSystem, constraint::Constraint)
	set_to_zero = system.K[constraint.index, :]
	for i in set_to_zero.nzind
		if i != constraint.index
			system.K[constraint.index, i] = 0.0
			system.K[i, constraint.index] = 0.0
		end
	end
	scale_by = system.K[constraint.index, constraint.index]
	system.K[constraint.index, constraint.index] *= constraint.coeff
	for (col_index, value) in constraint.constrain_to
		system.K[constraint.index, col_index] = -scale_by*value
		system.K[col_index, constraint.index] = -scale_by*value
	end
	system.F[constraint.index] = scale_by*constraint.inhomogeneity
end





# module boundary conditions ends here
end
# module boundary conditions ends here
