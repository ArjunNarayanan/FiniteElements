using Revise
using FiniteElements, isotropic, FileIO, LinearAlgebra

function body_force(x::AbstractArray)
    b = wave_number^2*(lambda + 2mu)*[sin(wave_number*x[1]),
                                       cos(wave_number*x[2])]
    return b
end

function displacement(x::AbstractArray)
	return [sin(wave_number*x[1]),
			cos(wave_number*x[2])]
end

function getDirichletData(mesh::Mesh, node_group::String)
	node_ids = mesh[:node_groups][node_group]
	dofs = [[1,2] for i in node_ids]
	vals = [[0.0,0.0] for i in node_ids]
	for (idx, id) in enumerate(node_ids)
		vals[idx] = displacement(mesh[:nodes][:,id])
	end
	return node_ids, dofs, vals
end

function getDirichletData(mesh::Mesh)
	groups = ["left", "bottom", "right", "top"]
	dirichlet_node_ids = Int[]
	dirichlet_node_dofs = Array{Array{Int64, 1}, 1}()
	dirichlet_values = Array{Array{Float64, 1}, 1}()
	for group in groups
		node_ids, dofs, vals = getDirichletData(mesh, group)
		append!(dirichlet_node_ids, node_ids)
		append!(dirichlet_node_dofs, dofs)
		append!(dirichlet_values, vals)
	end
	return dirichlet_node_ids, dirichlet_node_dofs, dirichlet_values
end

function applyDisplacementBCs(system::GlobalSystem, mesh::Mesh)
	dirichlet_node_ids, dirichlet_node_dofs, dirichlet_values =
			getDirichletData(mesh)
	applyDirichletBCs(system, dirichlet_node_ids, dirichlet_node_dofs,
		dirichlet_values)
end

function analyticalSolution(mesh::Mesh)
	nnodes = size(mesh[:nodes])[2]
	u = zeros(2, nnodes)
	for node_id in 1:nnodes
		u[:,node_id] = displacement(mesh[:nodes][:,node_id])
	end
	return u
end

const lambda = 10.0
const mu = 5.0
const wave_number = 2.0
const dofs = 2
q_order = 2
filename = "Quad_10.jld2"
folder = "test_elasticity/mesh/Quad/"
inputfile = folder*filename
mesh = load(inputfile)["mesh"]

system_matrix = SystemMatrix()
system_rhs = SystemRHS()
bilinearForm(λ -> lambda, μ -> mu, "body", mesh, q_order, system_matrix)
applyBodyForce(body_force, "body", mesh, q_order, system_rhs)
system = GlobalSystem(system_matrix, system_rhs, dofs)
applyDisplacementBCs(system, mesh)
solveDirect(system)

analytical_solution = analyticalSolution(mesh)
