using FileIO, Test
using FiniteElements, transformation_2d
using LinearAlgebra

##########################################################################
function readMesh(filename)
	mesh = load(filename)["mesh"]
	duplicateInterfaceNodes(mesh)
	haskey(mesh.data[:node_groups], "origin")
	haskey(mesh.data[:node_groups], "xmax")
	haskey(mesh.data[:node_groups], "xmax_core")
	haskey(mesh.data[:node_groups], "radial")
	haskey(mesh.data[:node_groups], "interface_core")
	haskey(mesh.data[:node_groups], "interface_shell")
	haskey(mesh.data[:node_groups], "shell_ir")
	haskey(mesh.data[:element_groups], "core")
	haskey(mesh.data[:element_groups], "shell")
	return mesh
end
##########################################################################

file = "Test.jld2"
mesh = readMesh(file)


##########################################################################
# Some tests
# First: Same number of elements in interface_Core and interface_shell
@test length(mesh.data[:element_groups]["interface_shell"][Line{2}]) ==
		length(mesh.data[:element_groups]["interface_core"][Line{2}])

# Second: All nodes in interface_core are in the core elements
I_core_elmt_ids = mesh.data[:element_groups]["interface_core"][Line{2}]
I_shell_elmt_ids = mesh.data[:element_groups]["interface_shell"][Line{2}]
core_elmt_ids = mesh.data[:element_groups]["core"][Triangle{3}]
shell_elmt_ids = mesh.data[:element_groups]["shell"][Triangle{3}]
I_core_elmts = mesh.data[:elements][Line{2}][:, I_core_elmt_ids]
I_shell_elmts = mesh.data[:elements][Line{2}][:, I_shell_elmt_ids]
core_elmts = mesh.data[:elements][Triangle{3}][:, core_elmt_ids]
shell_elmts = mesh.data[:elements][Triangle{3}][:, shell_elmt_ids]

for N in I_core_elmts
	@test N in core_elmts
	@test !(N in shell_elmts)
end

for N in I_shell_elmts
	@test N in shell_elmts
	@test !(N in core_elmts)
end


# Third: All nodes in I_core_elmts and I_shell elmts are at a 
# radius of 0.5
for i in eachindex(I_core_elmts)
	N1 = I_core_elmts[i]
	N2 = I_shell_elmts[i]
	node1 = mesh.data[:nodes][:, N1]
	node2 = mesh.data[:nodes][:, N2]
	@test abs(norm(node1) - 0.5) < 1e-10
	@test abs(norm(node2) - 0.5) < 1e-10
end









true
