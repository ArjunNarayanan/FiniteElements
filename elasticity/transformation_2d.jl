module transformation_2d

using FiniteElements, TensorOperations

export readMesh



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
	getFI(FI, ∇ϕI, E, ϵt)
Compute the contraction:
	∇sym(ϕI) * E * ϵt
This is the `[I]` entry of the element right-hand-side
vector. Here `∇sym` refers to the symmetrized gradient.
"""
function getFI(FI, ∇ϕI, E, ϵt)
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
						assembler, FI, E, ϵt)
	reinit(assembler)
	reinit(mapping, nodes)
	N = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		for I in 1:N
			∇ϕI = mapping[:gradients][I,q]
			getFI(FI, ∇ϕI, E, ϵt)
			assembler.element_rhs[I] += FI*mapping[:dx][q]
		end
	end
	updateSystem(assembler, node_ids)
end


"""
	assembleSystem(mesh::Mesh{spacedim}, Ec, Es, ϵts) where spacedim
Assemble the `GlobalSystem` for the current `mesh`. `Ec` and
`Es` are the 4th order elasticity tensors for the core and shell
respectively. `ϵts` is the transformation strain in the shell.
Currently only uses `Triangle{3}` elements, but can be generalized
in the future.
"""
function assembleSystem(mesh::Mesh{spacedim},
					Ec, Es, ϵts) where spacedim
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

	@time for elem_id in mesh.data[:element_groups]["core"]
		node_ids = mesh.data[:elements][elType][:,elem_id]
		nodes = mesh.data[:nodes][:, node_ids]
		assembleElementMatrix(node_ids, nodes, mapping,
					assembler, KIJ, Ec)
	end

	println("\tAssembling SHELL elements")

	@time for elem_id in mesh.data[:element_groups]["shell"]
		node_ids = mesh.data[:elements][elType][:,elem_id]
		nodes = mesh.data[:nodes][:, node_ids]
		assembleElementMatrix(node_ids, nodes, mapping,
					assembler, KIJ, Es)
		assembleElementRHS(node_ids, nodes, mapping,
					assembler, FI, Es, ϵts)
	end
	system = GlobalSystem(assembler)
	return system
end






# module ends here
end
# module ends here