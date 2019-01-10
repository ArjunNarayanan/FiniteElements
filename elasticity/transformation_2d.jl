module transformation_2d

using FiniteElements, TensorOperations, LinearAlgebra

export assembleSystem, duplicateInterfaceNodes

const spacedim = 2
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
	updateSystem(assembler, node_ids)
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
	updateSystem(assembler, node_ids)
end


"""
	assembleSystem(mesh::Mesh{spacedim}, Ec, Es, θ0) where spacedim
Assemble the `GlobalSystem` for the current `mesh`. `Ec` and
`Es` are the 4th order elasticity tensors for the core and shell
respectively. `θ0` is the volumetric transformation strain in the shell.
Currently only uses `Triangle{3}` elements, but can be generalized
in the future.
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

	elType = Triangle{3}
	dofs = 2

	total_ndofs = size(mesh.data[:nodes])[2]*dofs

	println("----------SOLVING FE PROBLEM----------")
	println("\tTotal Number of DOFs = ", total_ndofs)

	mapping = Map{elType,spacedim}(1, :gradients)
	assembler = Assembler(elType, dofs)

	KIJ = zeros(dofs,dofs)
	FI = zeros(dofs)

	println("\tAssembling CORE elements")

	@time for elem_id in mesh.data[:element_groups]["core"][elType]
		node_ids = mesh.data[:elements][elType][:,elem_id]
		nodes = mesh.data[:nodes][:, node_ids]
		assembleElementMatrix(node_ids, nodes, mapping,
					assembler, KIJ, Ec)
	end

	println("\tAssembling SHELL elements")

	@time for elem_id in mesh.data[:element_groups]["shell"][elType]
		node_ids = mesh.data[:elements][elType][:,elem_id]
		nodes = mesh.data[:nodes][:, node_ids]
		assembleElementMatrix(node_ids, nodes, mapping,
					assembler, KIJ, Es)
		assembleElementRHS(node_ids, nodes, mapping,
					assembler, FI, λs, μs, ϵt)
	end
	system = GlobalSystem(assembler)
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
		elmts = mesh.data[:elements][key][elm_ids]
		append!(node_ids, elmts[:])
	end
	unique!(node_ids)
	return node_ids
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
	n_1D_int_elmts_core = length(mesh.data[:element_groups]["interface_core"][Line{2}])
	##############################################

	##############################################
	shell_2D_elmt_ids = mesh.data[:element_groups]["shell"][Triangle{3}]
	shell_2D_elmts = mesh.data[:elements][Triangle{3}][:,shell_2D_elmt_ids]
	shell_1D_elmt_ids = range(n_1D_int_elmts_core+1, length = 2*n_1D_int_elmts_core)
	mesh.data[:element_groups]["interface_shell"] = Dict()
	mesh.data[:element_groups]["interface_shell"][Line{2}] = shell_1D_elmt_ids
	##############################################
end























# module ends here
end
# module ends here