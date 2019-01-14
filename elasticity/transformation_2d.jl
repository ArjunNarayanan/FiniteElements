module transformation_2d

using FiniteElements, TensorOperations, LinearAlgebra, SparseArrays

export assembleSystem, duplicateInterfaceNodes, 
		applyFreeSlipInterface

const spacedim = 2
const dofs = 2
δ = diagm(0 => ones(spacedim))


"""
	getKIJ(KIJ, ∇ϕI, E, ∇ϕJ)
Compute the contraction:
	∇sym(ϕI) * E * ∇sym(ϕJ) 
which is the `[I,J]` entry in the element stiffness 
matrix. Here `∇sym` refers to the symmetric gradient.
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
	assembleElementMatrix(nodes::Array{Float64, 2}, 
							   mapping::Map,
							   assembler::Assembler, 
							   KIJ::Array{Float64, 2}, 
							   E::Array{Float64, 4})
Compute the entries of the element matrix for the equilibrium equation
of linear elasticity. Store the entries into `assembler.element_matrix`.
"""
function assembleElementMatrix(nodes::Array{Float64, 2}, 
							   mapping::Map,
							   assembler::Assembler, 
							   KIJ::Array{Float64, 2}, 
							   E::Array{Float64, 4})

	reinit(assembler)
	reinit(mapping, nodes)

	Nnodes = length(mapping.master.basis.functions)

	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		for I in 1:Nnodes
			∇ϕI = mapping[:gradients][I,q]
			for J in 1:Nnodes
				∇ϕJ = mapping[:gradients][J,q]
				getKIJ(KIJ, ∇ϕI, E, ∇ϕJ)
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end


"""
	assembleElementRHS(nodes::Array{Float64, 2}, 
						    mapping::Map,
							assembler::Assembler, 
							FI::Array{Float64, 1}, λs::Float64, 
							μs::Float64, ϵt::Array{Float64, 2})
Compute the entries of the element RHS for the equilibrium equation
of linear elasticity driven by transformation strain `ϵt`. 
Store the entries into `assembler.element_rhs`.
"""
function assembleElementRHS(nodes::Array{Float64, 2}, 
						    mapping::Map,
							assembler::Assembler, 
							FI::Array{Float64, 1}, λs::Float64, 
							μs::Float64, ϵt::Array{Float64, 2})
	reinit(assembler)
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		for I in 1:Nnodes
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
function assembleSystem(mesh::Mesh{spacedim},
					λc, μc, λs, μs, θ0) where spacedim
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
			assembleElementMatrix(nodes, mapping,
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
			assembleElementMatrix(nodes, mapping,
						assembler, KIJ, Es)
			updateSystemMatrix(system_matrix, 
					assembler.element_matrix, node_ids,
					assembler.ndofs)
			assembleElementRHS(nodes, mapping,
						assembler, FI, λs, μs, ϵt)
			updateSystemRHS(system_rhs, assembler.element_rhs,
				node_ids, assembler.ndofs)
		end
	end

	system = GlobalSystem(system_matrix, system_rhs, dofs)
	return system
end

"""
	getNormal(T::Array{Float64, 1}, N::Array{Float64, 1})
Compute the outward normal `N` from the tangent `T`.
"""
function getNormal(T::Array{Float64, 1}, N::Array{Float64, 1})
	scale = sqrt(T[1]^2 + T[2]^2)
	N[1] = T[2]/scale
	N[2] = -T[1]/scale
end

"""
	getKIJ(KIJ::Array{Float64, 2}, N::Array{Float64, 1})
Compute the outer product `N` and store it in `KIJ`.
"""
function getKIJ(KIJ::Array{Float64, 2}, N::Array{Float64, 1})
	for i in 1:length(N)
		for j in 1:length(N)
			KIJ[i,j] = N[i]*N[j]
		end
	end
end

"""
	assembleElementMatrix(nodes::Array{Float64, 2}, 
						  mapping::Map, 
						  assembler::Assembler, 
						  KIJ::Array{Float64, 2},
						  KIJ_temp::Array{Float64, 2},
						  normal::Array{Float64, 1},
						  penalty::Float64)
Compute the entries of the element matrix for the free-slip interface
condition. Store the entries into `assembler.element_matrix`."""
function assembleElementMatrix(nodes::Array{Float64, 2}, 
							   mapping::Map, 
							   assembler::Assembler, 
							   KIJ::Array{Float64, 2},
							   KIJ_temp::Array{Float64, 2},
							   normal::Array{Float64, 1},
							   penalty::Float64)
	reinit(mapping, nodes)
	reinit(assembler)

	Nnodes = length(mapping.master.basis.functions)

	for q in eachindex(mapping.master.quadrature.points)
		getNormal(mapping[:jacobian][q][:,1], normal)
		getKIJ(KIJ, normal)
		KIJ *= penalty
		for I in 1:Nnodes
			ϕI = mapping[:values][I,q]
			for J in 1:Nnodes
				ϕJ = mapping[:values][J,q]
				KIJ_temp[:] = ϕI*ϕJ*KIJ[:]
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end

"""
	applyFreeSlipInterface(system::GlobalSystem, 
		mesh::Mesh; penalty = 1000)
Modifies `system` to ensure that degrees-of-freedom normal
to the interface have continuity via a penalty formulation.
"""
function applyFreeSlipInterface(system::GlobalSystem, 
			mesh::Mesh{spacedim}; penalty = 1000, q_order = 1) where spacedim
	ndofs = system.ndofs

	println("\tEnforcing free slip interface")
	elTypes = keys(mesh.data[:element_groups]["interface_core"])

	penalty = penalty*maximum(abs.(diag(system.K)))

	KIJ = zeros(ndofs, ndofs)
	KIJ_temp = zeros(ndofs, ndofs)
	normal = zeros(2)

	system_matrix = SystemMatrix()

	@time for elType in elTypes
		mapping = Map{elType,spacedim}(q_order, :values, :gradients)
		assembler = Assembler(elType, ndofs)
		for i in eachindex(mesh.data[:element_groups]["interface_core"][elType])
			c_el_id = mesh.data[:element_groups]["interface_core"][elType][i]
			s_el_id = mesh.data[:element_groups]["interface_shell"][elType][i]

			c_n_ids = mesh.data[:elements][elType][:, c_el_id]
			s_n_ids = mesh.data[:elements][elType][:, s_el_id]
			nodes = mesh.data[:nodes][:, c_n_ids]
			assembleElementMatrix(nodes, mapping, assembler, 
									KIJ, KIJ_temp, normal, penalty)
			updateSystemMatrix(system_matrix,
				assembler.element_matrix, c_n_ids,
				ndofs)
			updateSystemMatrix(system_matrix,
				assembler.element_matrix, s_n_ids, 
				ndofs)
			for j in 1:length(assembler.element_matrix)
				assembler.element_matrix[j] *= -1.0
			end
			updateSystemMatrix(system_matrix,
				assembler.element_matrix, s_n_ids,
				c_n_ids, ndofs)
			updateSystemMatrix(system_matrix,
				assembler.element_matrix, c_n_ids,
				s_n_ids, ndofs)
		end
	end
	n_g_dof = size(system.K)
	K2 = sparse(system_matrix.I, system_matrix.J,
				system_matrix.vals, n_g_dof[1], n_g_dof[2])
	system.K += K2
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
end






















# module ends here
end
# module ends here