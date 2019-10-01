using FileIO, WriteVTK, LinearAlgebra
using Revise
using FiniteElements, IsotropicElasticity

function body_force(x::AbstractArray)
    b = wave_number^2*(lambda + 2mu)*[sin(wave_number*x[1]),
                                       cos(wave_number*x[2])]
    return b
end

function analyticalDisplacement(x::AbstractArray)
	return [sin(wave_number*x[1]),
			cos(wave_number*x[2])]
end

function getDirichletData(mesh::Mesh, node_group::String)
	node_ids = mesh[:node_groups][node_group]
	dofs = [[1,2] for i in node_ids]
	vals = [[0.0,0.0] for i in node_ids]
	for (idx, id) in enumerate(node_ids)
		vals[idx] = analyticalDisplacement(mesh[:nodes][:,id])
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
		u[:,node_id] = analyticalDisplacement(mesh[:nodes][:,node_id])
	end
	return u
end

function solutionError(displacement::AbstractArray,
	analytical_solution::AbstractArray)

	err = 0.0
	for i in 1:size(displacement)[2]
		difference = displacement[:,i] - analytical_solution[:,i]
		err += norm(difference)/norm(analytical_solution[:,i])
	end
	return err
end

const lambda = 10.0
const mu = 5.0
const wave_number = 5.0
const dofs = 2
q_order = 2
filename = "Quad_100.jld2"
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

displacement = reshape(system.D, 2, :)
analytical_solution = analyticalSolution(mesh)
outfolder = "test_elasticity/Output/"
try
	mkpath(outfolder)
catch
	error("Failed to create output folder")
end

filename = split(filename, ".")[1]
outfilename = outfolder*filename
output = Output(outfilename, [Triangle{3}, Quadrilateral{4}], mesh)
writePointVectors(output, displacement, "displacement")
writePointVectors(output, analytical_solution, "analytical")
outfiles = vtk_save(output.vtkfile)

err = solutionError(displacement, analytical_solution)
println("Error = ", err)
