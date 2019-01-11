module transformation_2d

using FiniteElements, TensorOperations, LinearAlgebra

export assembleSystem, duplicateInterfaceNodes

const spacedim = 2
const dofs = 2
δ = diagm(0 => ones(spacedim))


"""
	getKIJ(KIJ, ∇ϕI, E, ∇ϕJ)
Compute the contraction:
	∇sym(ϕI) * E * ∇sym(ϕJ) 
which is the `[I,J]` entry in the element stiffness 
matrix. Here `∇sym` refers to the symmetrized gradient.
"""
function getKIJ(KIJ, ∇ϕI, E, ∇ϕJ)
	@tensor begin
		KIJ[p,r] = 0.5*(δ[i,p]*∇ϕI[j] + δ[j,p]*∇ϕI[i])*E[i,j,k,l]*0.5*(δ[k,r]*∇ϕJ[l] + δ[l,r]*∇ϕJ[k])
	end
end

"""
	getFI(FI, ∇ϕI, λs, μs, ϵt)
Compute the contraction:
	∇sym(ϕI) * E * ϵt
This is the `[I]` entry of the element right-hand-side
vector. Here `∇sym` refers to the symmetrized gradient.
"""
function getFI(FI, ∇ϕI, λs, μs, ϵt)
	θ0 = 3/2*tr(ϵt)
	@tensor begin
		FI[p] = 0.5*(δ[i,p]*∇ϕI[j] + δ[j,p]*∇ϕI[i])*(λs*θ0*δ[i,j] + 2μs*ϵt[i,j])
	end
end


"""
	assembleElementMatrix(node_ids, nodes, mapping,
							assembler, KIJ, E, ϵt)
Compute the entries of the element matrix. Assemble
the entries into `assembler.system_matrix`.
"""
function assembleElementMatrix(node_ids, nodes, mapping,
							assembler, KIJ, E)

	reinit(assembler)
	reinit(mapping, nodes)

	N = length(mapping.master.basis.functions)
	# Loop over quadrature points
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		for I in 1:N
			∇ϕI = mapping[:gradients][I,q]
			for J in 1:N
				∇ϕJ = mapping[:gradients][J,q]
				getKIJ(KIJ, ∇ϕI, E, ∇ϕJ)
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end


"""
	assembleElementRHS(node_ids, nodes, mapping,
						assembler, FI, E, ϵt)
Compute the entries of the element RHS. Assemble the
entries into `assembler.system_rhs`.
"""
function assembleElementRHS(node_ids, nodes, mapping,
						assembler, FI, λs, μs, ϵt)
	reinit(assembler)
	reinit(mapping, nodes)
	N = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		for I in 1:N
			∇ϕI = mapping[:gradients][I,q]
			getFI(FI, ∇ϕI, λs, μs, ϵt)
			assembler.element_rhs[I] += FI*mapping[:dx][q]
		end
	end
end


"""
	assembleSystem(mesh::Mesh{spacedim}, Ec, Es, θ0) where spacedim
Assemble the `GlobalSystem` for the current `mesh`. `Ec` and
`Es` are the 4th order elasticity tensors for the core and shell
respectively. `θ0` is the volumetric transformation strain in the shell.
"""
function assembleSystem(mesh::Mesh{2},
					λc, μc, λs, μs, θ0)
	Ec = zeros(ntuple(x -> spacedim, 4)...)
	Es = zeros(ntuple(x -> spacedim, 4)...)

	ϵt = θ0/3*δ

	@tensor begin
		Ec[i,j,k,l] = λc*δ[i,j]*δ[k,l] + μc*(δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k])
	end
	@tensor begin
		Es[i,j,k,l] = λs*δ[i,j]*δ[k,l] + μs*(δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k])
	end

	core_elTypes = keys(mesh.data[:element_groups]["core"])
	shell_elTypes = keys(mesh.data[:element_groups]["shell"])
	elTypes = union(core_elTypes, shell_elTypes)

	total_ndofs = size(mesh.data[:nodes])[2]*dofs

	println("----------SOLVING FE PROBLEM----------")
	println("\tTotal Number of DOFs = ", total_ndofs)

	mapping_dict = Dict()
	assembler_dict = Dict()
	
	for elType in elTypes
		mapping_dict[elType] = Map{elType,spacedim}(1, :gradients)
		assembler_dict[elType] = Assembler(elType, dofs)
	end

	KIJ = zeros(dofs,dofs)
	FI = zeros(dofs)

	system_matrix = SystemMatrix()
	system_rhs = SystemRHS()

	println("\tAssembling CORE elements")

	@time for elType in core_elTypes
		mapping = mapping_dict[elType]
		assembler = assembler_dict[elType]
		for elem_id in mesh.data[:element_groups]["core"][elType]
			node_ids = mesh.data[:elements][elType][:,elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			assembleElementMatrix(node_ids, nodes, mapping,
						assembler, KIJ, Ec)
			updateSystemMatrix(system_matrix, 
				assembler.element_matrix, node_ids, 
				assembler.ndofs)
		end
	end

	println("\tAssembling SHELL elements")

	@time for elType in shell_elTypes
		mapping = mapping_dict[elType]
		assembler = assembler_dict[elType]
		for elem_id in mesh.data[:element_groups]["shell"][elType]
			node_ids = mesh.data[:elements][elType][:,elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			assembleElementMatrix(node_ids, nodes, mapping,
						assembler, KIJ, Es)
			updateSystemMatrix(system_matrix, 
					assembler.element_matrix, node_ids,
					assembler.ndofs)
			assembleElementRHS(node_ids, nodes, mapping,
						assembler, FI, λs, μs, ϵt)
			updateSystemRHS(system_rhs, assembler.element_rhs,
				node_ids, assembler.ndofs)
		end
	end

	system = GlobalSystem(system_matrix, system_rhs, dofs)
	return system
end

"""
	getInterfaceNodeIDs(mesh)
Return a list of all unique node IDs of nodes on the interface.
"""
function getInterfaceNodeIDs(mesh)
	node_ids = Int[]
	for key in keys(mesh.data[:element_groups]["interface"])
		elm_ids = mesh.data[:element_groups]["interface"][key]
		elmts = mesh.data[:elements][key][:, elm_ids]
		append!(node_ids, elmts[:])
	end
	unique!(node_ids)
	return node_ids
end

"""
	replaceNodeIDs(input_array, old_node_ids, new_node_ids)
Replace all occurrences of node IDs in `old_node_ids` with their 
corresponding values in `new_node_ids`.
"""
function replaceNodeIDs!(input_array, old_node_ids, new_node_ids)
	@assert length(old_node_ids) == length(new_node_ids)
	for i in eachindex(old_node_ids)
		N_old = old_node_ids[i]
		N_new = new_node_ids[i]
		replace!(input_array, N_old => N_new)
	end
end


"""
	duplicateInterfaceNodes(mesh::Mesh{2})
Duplicate the nodes on the interface. Modify the connectivity
of shell elements to refer to these duplicated nodes. Used to 
model an incoherent interface.
"""
function duplicateInterfaceNodes(mesh::Mesh{2})
	##############################################
	# First: duplicate all the nodes on the interface
	n_old_nodes = size(mesh.data[:nodes])[2]
	old_int_node_ids = getInterfaceNodeIDs(mesh)
	interface_nodes = mesh.data[:nodes][:,old_int_node_ids]
	mesh.data[:nodes] = hcat(mesh.data[:nodes], interface_nodes)
	n_add_nodes = length(old_int_node_ids)
	new_int_node_ids = range(n_old_nodes+1, length = n_add_nodes)
	##############################################

	##############################################
	# Second: redefine mesh.data[:element_groups]["interface"]
	mesh.data[:element_groups]["interface_core"] = pop!(mesh.data[:element_groups], "interface")
	##############################################

	##############################################
	# Third: Replace all occurrences of the old node ids in the shell
	# elements with the new node ids
	for key in keys(mesh.data[:element_groups]["shell"])
		shell_2D_elmt_ids = mesh.data[:element_groups]["shell"][key]
		shell_2D_elmts = mesh.data[:elements][key][:,shell_2D_elmt_ids]
		replaceNodeIDs!(shell_2D_elmts, old_int_node_ids, new_int_node_ids)
		mesh.data[:elements][key][:,shell_2D_elmt_ids] = shell_2D_elmts
	end
	##############################################

	##############################################
	# Fourth: add the additional 1D elements for the new nodes 
	# into mesh.data[:elements]
	mesh.data[:element_groups]["interface_shell"] = Dict()
	for key in keys(mesh.data[:element_groups]["interface_core"])
		core_1D_elmt_ids = mesh.data[:element_groups]["interface_core"][key]
		core_1D_elmts = mesh.data[:elements][key][:, core_1D_elmt_ids]
		shell_1D_elmts = copy(core_1D_elmts)
		replaceNodeIDs!(shell_1D_elmts, old_int_node_ids, new_int_node_ids)
		e_last = size(mesh.data[:elements][key])[2]
		mesh.data[:elements][key] = hcat(mesh.data[:elements][key], shell_1D_elmts)
		n_shell_1D_elmts = size(shell_1D_elmts)[2]
		shell_1D_elmt_ids = collect(range(e_last+1, length = n_shell_1D_elmts))
		mesh.data[:element_groups]["interface_shell"][key] = shell_1D_elmt_ids
	end
	##############################################

	##############################################
	# Fifth: Find the node corresponding to xmax_core
	# get the shell equivalent of this node, and store it in
	# shell_ir
	xmax_core_id = mesh.data[:node_groups]["xmax_core"][1]
	id = findfirst(x -> x == xmax_core_id, old_int_node_ids)
	shell_ir_id = new_int_node_ids[id]
	mesh.data[:node_groups]["shell_ir"] = [shell_ir_id]
	##############################################	
end






















# module ends here
end
# module ends here