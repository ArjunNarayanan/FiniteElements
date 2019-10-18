module transformation_2d

using FiniteElements, TensorOperations, LinearAlgebra, SparseArrays

export assembleSystem, duplicateInterfaceNodes,
		applyFreeSlipInterface, computeElementAveragedStrain,
		computeElementAveragedStress, writeCellStrainComponents,
		writeCellStressComponents

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

function getKIJ(KIJ, ∇ϕI, ∇ϕJ, λ, μ)
	@tensor begin
		KIJ[p,r] = 0.5*(δ[i,p]*∇ϕI[j] + δ[j,p]*∇ϕI[i])*
						(λ*δ[i,j]*δ[k,l] + 2*μ*δ[i,k]*δ[j,l])*
						0.5*(δ[k,r]*∇ϕJ[l] + δ[l,r]*∇ϕJ[k])
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
	function assembleElementMatrix(nodes::Array{Float64, 2},
		mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2},
		E::Array{Float64, 4})
Compute the entries of the element matrix for the equilibrium equation
of linear elasticity. Store the entries into `assembler.element_matrix`.
"""
function assembleElementMatrix(nodes::Array{Float64, 2},
	mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2},
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
	function assembleElementMatrix(nodes::Array{Float64, 2},
		mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2}, λ::Function,
		μ::Function)
Compute the entries of the element matrix for the equilibrium equation of
linear elasticity. Lame coefficients `λ(X), μ(X)` are computed at the
coordinates `X`.
"""
function assembleElementMatrix(nodes::Array{Float64, 2},
	mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2}, λ::Function,
	μ::Function)

	reinit(assembler)
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)

	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		Xq = mapping[:coordinates][q]
		λq = λ(Xq)
		μq = μ(Xq)
		for I in 1:Nnodes
			∇ϕI = mapping[:gradients][I,q]
			for J in 1:Nnodes
				∇ϕJ = mapping[:gradients][J,q]
				getKIJ(KIJ, ∇ϕI, ∇ϕJ, λq, μq)
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end


"""
	function assembleElementRHS(nodes::Array{Float64, 2}, mapping::Map,
		assembler::Assembler, FI::Array{Float64, 1}, λs::Float64, μs::Float64,
		ϵt::Array{Float64, 2})
Compute the entries of the element RHS for the equilibrium equation
of linear elasticity driven by transformation strain `ϵt`.
Store the entries into `assembler.element_rhs`.
"""
function assembleElementRHS(nodes::Array{Float64, 2}, mapping::Map,
	assembler::Assembler, FI::Array{Float64, 1}, λs::Float64, μs::Float64,
	ϵt::Array{Float64, 2})

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
	function assembleElementRHS(nodes::Array{Float64, 2}, mapping::Map,
		asembler::Assembler, FI::Array{Float64, 1}, λ::Function, μ::Function,
		ϵt::Array{Float64, 2})
Compute the entries of the element RHS for the equilibrium equation of linear
elasticity driven by transformation strain `ϵt`. The Lame coefficients
`λ(X),μ(X)` are computed at quadrature point `X`.
"""
function assembleElementRHS(nodes::Array{Float64, 2}, mapping::Map,
	assembler::Assembler, FI::Array{Float64, 1}, λ::Function, μ::Function,
	ϵt::Array{Float64, 2})

	reinit(assembler)
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		Xq = mapping[:coordinates][q]
		λq = λ(Xq)
		μq = μ(Xq)
		for I in 1:Nnodes
			∇ϕI = mapping[:gradients][I,q]
			getFI(FI, ∇ϕI, λq, μq, ϵt)
			assembler.element_rhs[I] += FI*mapping[:dx][q]
		end
	end
end


"""
	function assembleSystem(mesh::Mesh{spacedim},
		λc::Float64, μc::Float64, λs::Float64, μs::Float64, θ0::Float64;
		q_order = 1) where spacedim
Assemble the `GlobalSystem` for the current `mesh`. `λc, μc` are the Lame coefficients
in the core, and `λs, μs` are the Lame coefficients in the shell.
`θ0` is the volumetric transformation strain in the shell. `q_order` is the quadrature
order to be used.
"""
function assembleSystem(mesh::Mesh{spacedim},
	λc::Float64, μc::Float64, λs::Float64, μs::Float64, θ0::Float64;
	q_order = 1) where spacedim

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
		mapping_dict[elType] = Map{elType,spacedim}(q_order, :gradients)
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
			node_ids = mesh.data[:elements][elType][:, elem_id]
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
			node_ids = mesh.data[:elements][elType][:, elem_id]
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
	function assembleSystem(mesh::Mesh{spacedim}, λ::Function, μ::Function,
		θ0::Float64; q_order = 1) where spacedim
Assemble the `GlobalSystem` for the current `mesh`. The Lame coefficients
`λ, μ` are given as functions over the domain. `θ0` is the volumetric
transformation strain in the shell. `q_order` is the quadrature order to be
used.
"""
function assembleSystem(mesh::Mesh{spacedim}, λ::Function, μ::Function,
	θ0::Float64; q_order = 1) where spacedim

	ϵt = θ0/3*δ

	core_elTypes = keys(mesh.data[:element_groups]["core"])
	shell_elTypes = keys(mesh.data[:element_groups]["shell"])
	elTypes = union(core_elTypes, shell_elTypes)

	total_ndofs = size(mesh.data[:nodes])[2]*dofs

	println("----------SOLVING FE PROBLEM----------")
	println("\tTotal Number of DOFs = ", total_ndofs)

	mapping_dict = Dict()
	assembler_dict = Dict()

	for elType in elTypes
		mapping_dict[elType] = Map{elType,spacedim}(q_order, :coordinates,
			:gradients)
		assembler_dict[elType] = Assembler(elType, dofs)
	end

	KIJ = zeros(dofs, dofs)
	FI = zeros(dofs)

	system_matrix = SystemMatrix()
	system_rhs = SystemRHS()

	println("\tAssembling CORE elements")

	@time for elType in core_elTypes
		mapping = mapping_dict[elType]
		assembler = assembler_dict[elType]
		for elem_id in mesh.data[:element_groups]["core"][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			assembleElementMatrix(nodes, mapping, assembler, KIJ, λ, μ)
			updateSystemMatrix(system_matrix, assembler.element_matrix,
				node_ids, assembler.ndofs)
		end
	end

	println("\tAssembling SHELL elements")

	@time for elType in shell_elTypes
		mapping = mapping_dict[elType]
		assembler = assembler_dict[elType]
		for elem_id in mesh.data[:element_groups]["shell"][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			assembleElementMatrix(nodes, mapping, assembler, KIJ, λ, μ)
			updateSystemMatrix(system_matrix, assembler.element_matrix,
				node_ids, assembler.ndofs)
			assembleElementRHS(nodes, mapping, assembler, FI, λ, μ, ϵt)
			updateSystemRHS(system_rhs, assembler.element_rhs, node_ids,
				assembler.ndofs)
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
				assembler.element_matrix[I,J] += KIJ_temp*mapping[:dx][q]
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

	println("\tEnforcing FREE SLIP interface condition")
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

"""
	updateStrain(strain::Array{Float64, 2},
			     elem_id::Int64,
			     node_ids::Array{Int64, 1},
			     nodes::Array{Float64, 2},
			     mapping::Map,
			     displacement::Array{Float64, 2})
Update `strain[:,elem_id]` with the corresponding
strain components.
"""
function updateStrain(strain::Array{Float64, 2},
					  elem_id::Int64,
					  node_ids::Array{Int64, 1},
					  nodes::Array{Float64, 2},
					  mapping::Map,
					  displacement::Array{Float64, 2})
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)

	for q in eachindex(mapping.master.quadrature.points)
		w = mapping.master.quadrature.weights[q]/sum(mapping.master.quadrature.weights)
		for I in 1:Nnodes
			strain[1, elem_id] += mapping[:gradients][I,q][1]*displacement[1,node_ids[I]]*w
			strain[2, elem_id] += 0.5*( mapping[:gradients][I,q][2]*displacement[1,node_ids[I]] +
										mapping[:gradients][I,q][q]*displacement[2,node_ids[I]] )*w
			strain[3, elem_id] += mapping[:gradients][I,q][2]*displacement[2,node_ids[I]]*w
		end
	end
end

"""
	computeElementAveragedStrain(mesh::Mesh{spacedim},
									  displacement::Array{Float64, 2})
Returns the element averaged strain in a dictionary whose keys are element types,
and the associated values are arrays of size `(3, n_elements)` with entries:
	[e11 ...
	 e12 ...
	 e22 ...]
`displacement` is expected to be a `(2,n_nodes)` dimension array with displacement
components.
"""
function computeElementAveragedStrain(mesh::Mesh,
									  displacement::Array{Float64, 2};
									  q_order = 1)
	core_elTypes = keys(mesh.data[:element_groups]["core"])
	shell_elTypes = keys(mesh.data[:element_groups]["shell"])
	elTypes = union(core_elTypes, shell_elTypes)

	mapping_dict = Dict()
	strain_dict = Dict()

	for elType in elTypes
		mapping_dict[elType] = Map{elType,spacedim}(q_order, :gradients)
		n_elements = size(mesh.data[:elements][elType])[2]
		strain_dict[elType] = zeros(3, n_elements)
	end

	println("\tComputing element-averaged strain for CORE")

	@time for elType in core_elTypes
		mapping = mapping_dict[elType]
		strain = strain_dict[elType]
		for elem_id in mesh.data[:element_groups]["core"][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			updateStrain(strain, elem_id, node_ids, nodes, mapping, displacement)
		end
	end
	println("\tComputing element-averaged strain for SHELL")
	@time for elType in shell_elTypes
		mapping = mapping_dict[elType]
		strain = strain_dict[elType]
		for elem_id in mesh.data[:element_groups]["shell"][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			updateStrain(strain, elem_id, node_ids, nodes, mapping, displacement)
		end
	end
	return strain_dict
end


function e11(strain::Array{Float64, 2})
	return strain[1,:]
end

function e12(strain::Array{Float64, 2})
	return strain[2,:]
end

function e22(strain::Array{Float64, 2})
	return strain[3,:]
end

function volumetric(strain::Array{Float64, 2})
	return e11(strain) + e22(strain)
end

function dev_e11(strain::Array{Float64, 2})
	data = e11(strain) - 1/3*volumetric(strain)
	return data
end

function dev_e22(strain::Array{Float64, 2})
	data = e22(strain) - 1/3*volumetric(strain)
	return data
end

"""
	writeCellStrainComponents(output::Output, strain_dict::Array{Float64, 2},
									args::Vararg{Symbol})
Write the strain components specified in `args` into `output.vtkfile`. The data
is treated as cell data. See `computeElementAveragedStrain` to compute `strain_dict`.
# Supported arguments for `args`:
- `:e11`
- `:e12`
- `:e22`
- `:volumetric`
- `:dev_e11` - [1,1] component of deviatoric strain
- `:dev_e22` - [2,2] component of deviatoric strain
"""
function writeCellStrainComponents(output::Output, strain_dict::Dict,
									args::Vararg{Symbol})
	@assert !isempty(args) "Argument list cannot be empty"
	strain = hcat([strain_dict[e] for e in output.elTypes]...)
	for arg in args
		data = eval(arg)(strain)
		writeCellScalars(output, data, string(arg))
	end
end


function s11(stress::Array{Float64, 2})
	return stress[1,:]
end

function s12(stress::Array{Float64, 2})
	return stress[2,:]
end

function s22(stress::Array{Float64, 2})
	return stress[3,:]
end

function s33(stress::Array{Float64, 2})
	return stress[4,:]
end

function pressure(stress::Array{Float64, 2})
	return -1/3*(s11(stress) + s22(stress) + s33(stress))
end

function dev_s11(stress::Array{Float64, 2})
	return s11(stress) + pressure(stress)
end

function dev_s22(stress::Array{Float64, 2})
	return s22(stress) + pressure(stress)
end

function dev_s33(stress::Array{Float64, 2})
	return s33(stress) + pressure(stress)
end

function dev_s_norm(stress::Array{Float64, 2})
	s_norm = sqrt.( dev_s11(stress).^2 +
					2*s12(stress).^2 +
					dev_s22(stress).^2 +
					dev_s33(stress).^2 )
	return s_norm
end


"""
	writeCellStressComponents(output::Output, stress_dict::Dict,
									args::Vararg{Symbol})
Write the stress components specified in `args` into `output.vtkfile`. The data
is treated as cell data. See `computeElementAveragedStress` to compute `stress_dict`.
# Supported arguments for `args`:
- `:s11`
- `:s12`
- `:s22`
- `:s33`
- `:pressure`
- `:dev_s11` - [1,1] component of deviatoric stress
- `:dev_s22` - [2,2] component of deviatoric stress
- `:dev_s33` - [3,3] component of deviatoric stress
- `:dev_s_norm` - norm of the deviatoric stress
"""
function writeCellStressComponents(output::Output, stress_dict::Dict,
									args::Vararg{Symbol})
	@assert !isempty(args) "Argument list cannot be empty"
	stress = hcat([stress_dict[e] for e in output.elTypes]...)
	for arg in args
		data = eval(arg)(stress)
		writeCellScalars(output, data, string(arg))
	end
end




"""
	updateStress(stress::SubArray, strain::SubArray, λ::Float64, μ::Float64,
					θ0::Float64)
Update `stress` using a plane strain linear elastic constitutive model with `λ` and `μ` as
the usual Lame parameters. `θ0` is taken as the stress free volumetric transformation
strain. The strain components are assumed as follows:
	strain[1] = ϵ[1,1]
	strain[2] = ϵ[1,2]
	strain[3] = ϵ[2,2]
The stress components are:
	stress[1] = σ[1,1]
	stress[2] = σ[1,2]
	stress[3] = σ[2,2]
	stress[4] = σ[3,3]
"""
function updateStress(stress::SubArray, strain::SubArray, λ::Float64, μ::Float64,
						θ0::Float64)
	stress[1] = λ*(strain[1] + strain[3] - θ0) + 2μ*(strain[1] - θ0/3)
	stress[2] = 2μ*strain[2]
	stress[3] = λ*(strain[1] + strain[3] - θ0) + 2μ*(strain[3] - θ0/3)
	stress[4] = λ*(strain[1] + strain[3] - θ0) + 2μ*(-θ0/3)
end


"""
	computeElementAveragedStress(mesh::Mesh,
			strain_dict::Dict, λc::Float64,
			μc::Float64, λs::Float64, μs::Float64, θ0::Array{Float64, 2})
Returns the element averaged stress in a dictionary whose keys are element types.
The associated values are arrays of size `(4,n_elements)`, where `n_elements` is the
number of elements of that type, with entries:
	[s11 ...
	 s12 ...
	 s22 ...
	 s33 ...]
`strain_dict` is expected to be a dict whose keys are element types. The value of each
key is an array of size `(3,n_elements)` where `n_elements` is the number of elements
of that type. The components of the array are expected to be:
	[e11 ...
	 e12 ...
	 e22 ...]
"""
function computeElementAveragedStress(mesh::Mesh,
			strain_dict::Dict, λc::Float64,
			μc::Float64, λs::Float64, μs::Float64, θ0::Float64)
	stress_dict = Dict()

	for etype in keys(strain_dict)
		n_elements = size(strain_dict[etype])[2]
		stress_dict[etype] = zeros(4, n_elements)
	end

	println("Computing element averaged stress for CORE")

	@time for etype in keys(mesh.data[:element_groups]["core"])
		strain = strain_dict[etype]
		stress = stress_dict[etype]
		for elem_id in mesh.data[:element_groups]["core"][etype]
			sigma = view(stress, :, elem_id)
			epsilon = view(strain, :, elem_id)
			updateStress(sigma, epsilon, λc, μc, 0.0)
		end
	end

	println("Computing element averaged stress for SHELL")

	@time for etype in keys(mesh.data[:element_groups]["shell"])
		strain = strain_dict[etype]
		stress = stress_dict[etype]
		for elem_id in mesh.data[:element_groups]["shell"][etype]
			sigma = view(stress, :, elem_id)
			epsilon = view(strain, :, elem_id)
			updateStress(sigma, epsilon, λs, μs, θ0)
		end
	end

	return stress_dict
end










# module ends here
end
# module ends here
