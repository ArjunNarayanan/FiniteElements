module isotropic

using FiniteElements, TensorOperations

"""
	getKIJ(KIJ::Array{Float64, 2}, ∇ϕI::Array{Float64, 1},
		∇ϕJ::Array{Float64, 1}, λ::Float64, μ::Float64, δ::Array{Float64, 2})
Compute the weak form of linear elasticity at a single quadrature point. The
result is stored in `KIJ`.
"""
function getKIJ(KIJ::Array{Float64, 2}, ∇ϕI::Array{Float64, 1},
	∇ϕJ::Array{Float64, 1}, λ::Float64, μ::Float64, δ::Array{Float64, 2})

	@tensor begin
		KIJ[p,r] = 0.5*(δ[i,p]*∇ϕI[j] + δ[j,p]*∇ϕI[i])*
						(λ*δ[i,j]*δ[k,l] + 2*μ*δ[i,k]*δ[j,l])*
						0.5*(δ[k,r]*∇ϕJ[l] + δ[l,r]*∇ϕJ[k])
	end
end

"""
	function assembleElementMatrix(nodes::Array{Float64, 2},
		mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2},
		lambda::Function, mu::Function)
Compute the entries of the element matrix for the equilibrium equation of
linear elasticity. Lame coefficients `lambda(Xq), mu(Xq)` are computed at the
gauss quadrature points `Xq`.
"""
function assembleElementMatrix(nodes::Array{Float64, 2},
	mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2},
	lambda::Function, mu::Function, kronecker_delta::Array{Float64, 2})

	reinit(assembler)
	reinit(mapping, nodes)

	Nnodes = length(mapping.master.basis.functions)

	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		Xq = mapping[:coordinates][q]
		λq = lambda(Xq)
		μq = mu(Xq)
		for I in 1:Nnodes
			∇ϕI = mapping[:gradients][I,q]
			for J in 1:Nnodes
				∇ϕJ = mapping[:gradients][J,q]
				getKIJ(KIJ, ∇ϕI, ∇ϕJ, λq, μq, kronecker_delta)
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end


"""
    bilinearForm(element_group::String, lambda::Function, mu::Function,
        mesh::Mesh{spacedim}, system_matrix::SystemMatrix; q_order = 1) where spacedim
assemble the bilinear form for isotropic linear elasticity on the domain
specified by `element_group`. `lambda` and `mu` are the Lame coefficients
for an isotropic material.
"""
function bilinearForm(element_group::String, lambda::Function, mu::Function,
    mesh::Mesh{spacedim}, system_matrix::SystemMatrix; q_order = 1) where spacedim

	kronecker_delta = diagm(0 => ones(spacedim))
    elTypes = keys(mesh.data[:element_groups][element_group])

    KIJ = zeros(spacedim,spacedim)
    FI = zeros(spacedim)

    for elType in elTypes
        mapping = Map{elType,spacedim}(q_order, :gradients)
        assembler = Assembler(elType, spacedim)

        for elem_id in mesh.data[:element_groups][element_group][elType]
            node_ids = mesh.data[:elements][elType][:, elem_id]
            nodes = mesh.data[:nodes][:, node_ids]
            assembleElementMatrix(nodes, mapping, assembler, KIJ, lambda, mu,
				kronecker_delta)
			updateSystemMatrix(system_matrix, assembler.element_matrix,
				node_ids, assembler.ndofs)
        end
    end
end



function assembleElementRHS(nodes::Array{Float64, 2}, mapping::Map,
	assembler::Assembler, FI::Array{Float64, 1}, traction::Function)

	reinit(assembler)
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		Xq = mapping[:coordinates][q]
		tau = traction(Xq)
		for I in 1:Nnodes
			ϕI = mapping[:values][I,q]
			FI[:] = ϕI*tau[:]
			assembler.element_rhs[I] += FI*mapping[:dx][q]
		end
	end
end

"""
	linearForm(element_group::String, traction::Function,
		mesh::Mesh{spacedim}, system_rhs::SystemRHS; q_order = 1) where spacedim
assemble the linear form (i.e. right-hand-side) of linear elasticity.
`traction` is the load to be applied on the domain specified by `element_group`.
"""
function linearForm(element_group::String, traction::Function,
	mesh::Mesh{spacedim}, system_rhs::SystemRHS; q_order = 1) where spacedim

	elTypes = keys(mesh.data[:element_groups][element_group])
	FI = zeros(spacedim)

	for elType in elTypes
		mapping = Map{elType,spacedim}(q_order, :values, :gradients)
		assembler = Assembler(elType, spacedim)
		for elem_id in mesh.data[:element_groups][element_group][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			assembleElementRHS(nodes, mapping, assembler, FI, traction)
			updateSystemRHS(system_rhs, assembler.element_rhs, node_ids,
				assembler.ndofs)
		end
	end
end


# module isotropic ends here
end
# module isotropic ends here
