module isotropic

using FiniteElements, TensorOperations

"""
	bilinearForm(KIJ::Array{Float64, 2}, ∇ϕI::Array{Float64, 1},
		∇ϕJ::Array{Float64, 1}, λ::Float64, μ::Float64, δ::Array{Float64, 2})
Compute the bilinear form of linear elasticity at a single quadrature point.
The result is stored in `KIJ`.
"""
function bilinearForm(KIJ::Array{Float64, 2}, ∇ϕI::Array{Float64, 1},
	∇ϕJ::Array{Float64, 1}, λ::Float64, μ::Float64, δ::Array{Float64, 2})

	@tensor begin
		KIJ[p,r] = 0.5*(δ[i,p]*∇ϕI[j] + δ[j,p]*∇ϕI[i])*
						(λ*δ[i,j]*δ[k,l] + 2*μ*δ[i,k]*δ[j,l])*
						0.5*(δ[k,r]*∇ϕJ[l] + δ[l,r]*∇ϕJ[k])
	end
end

"""
	function bilinearForm(nodes::Array{Float64, 2},
		mapping::Map, assembler::Assembler, KIJ::Array{Float64, 2},
		lambda::Function, mu::Function)
Compute the bilinear form for the equilibrium equation of linear elasticity
on a particular element. Lame coefficients `lambda(Xq), mu(Xq)` are computed
at the gauss quadrature points `Xq`.
"""
function bilinearForm(nodes::Array{Float64, 2},
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
				bilinearForm(KIJ, ∇ϕI, ∇ϕJ, λq, μq, kronecker_delta)
				assembler.element_matrix[I,J] += KIJ*mapping[:dx][q]
			end
		end
	end
end


"""
    bilinearForm(element_group::String, lambda::Function, mu::Function,
        mesh::Mesh{spacedim}, system_matrix::SystemMatrix; q_order = 1) where spacedim
Assemble the bilinear form for isotropic linear elasticity on the domain
specified by `element_group`. `lambda` and `mu` are the Lame coefficients
for an isotropic material. The result is updated into the `system_matrix`.
"""
function bilinearForm(lambda::Function, mu::Function, element_group::String,
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
            bilinearForm(nodes, mapping, assembler, KIJ, lambda, mu,
				kronecker_delta)
			updateSystemMatrix(system_matrix, assembler.element_matrix,
				node_ids, assembler.ndofs)
        end
    end
end


"""
	linearForm(force::Function, nodes::Array{Float64, 2}, mapping::Map,
		assembler::Assembler, FI::Array{Float64, 1})
Assemble the right-hand-side vector of linear elasticity on a particular
element specified by `nodes`, `mapping`, and `Assembler`.
`force` is evaluated at each gauss quadrature point to obtain the force vector.
The result is stored in `FI`.
"""
function linearForm(force::Function, nodes::Array{Float64, 2}, mapping::Map,
	assembler::Assembler, FI::Array{Float64, 1})

	reinit(assembler)
	reinit(mapping, nodes)
	Nnodes = length(mapping.master.basis.functions)
	for q in eachindex(mapping.master.quadrature.points)
		(pq, wq) = mapping.master.quadrature[q]
		Xq = mapping[:coordinates][q]
		tau = force(Xq)
		for I in 1:Nnodes
			ϕI = mapping[:values][I,q]
			FI[:] = ϕI*tau[:]
			assembler.element_rhs[I] += FI*mapping[:dx][q]
		end
	end
end

"""
	linearForm(force::Function, element_group::String,
		mesh::Mesh{spacedim}, q_order::Int64,
		system_rhs::SystemRHS) where spacedim
assemble the linear form (i.e. right-hand-side) of linear elasticity.
`force` is the load to be applied on the domain specified by `element_group`.
`q_order` specifies the quadrature order to be used. The result is updated into
`system_rhs`.
"""
function linearForm(force::Function, element_group::String,
	mesh::Mesh{spacedim}, q_order::Int64, system_rhs::SystemRHS) where spacedim

	elTypes = keys(mesh.data[:element_groups][element_group])
	FI = zeros(spacedim)

	for elType in elTypes
		mapping = Map{elType,spacedim}(q_order, :values, :gradients)
		assembler = Assembler(elType, spacedim)
		for elem_id in mesh.data[:element_groups][element_group][elType]
			node_ids = mesh.data[:elements][elType][:, elem_id]
			nodes = mesh.data[:nodes][:, node_ids]
			linearForm(traction, nodes, mapping, assembler, FI)
			updateSystemRHS(system_rhs, assembler.element_rhs, node_ids,
				assembler.ndofs)
		end
	end
end

"""
	applyTraction(traction::Function, element_group::String,
		mesh::Mesh{spacedim}, q_order::Int64, system_rhs::SystemRHS)
Treats the function `traction` as a traction force applied on the group
specified by `element_group`. This is just a convenience function to
distinguish from a "body force" type of load.
"""
function applyTraction(traction::Function, element_group::String,
	mesh::Mesh{spacedim}, q_order::Int64, system_rhs::SystemRHS)

	linearForm(traction, element_group, mesh, q_order, system_rhs)
end

"""
	applyBodyForce(body_force::Function, element_group::String,
		mesh::Mesh{spacedim}, q_order::Int64, system_rhs::SystemRHS)
Treats the function `body_force` as a body force applied on the group
specified by `element_group`. In linear elasticity, the body force takes on a
negative sign when applied to the right-hand-side vector of the linear system.
"""
function applyBodyForce(body_force::Function, element_group::String,
	mesh::Mesh{spacedim}, q_order::Int64, system_rhs::SystemRHS)

	linearForm(x -> -body_force(x), element_group, mesh, q_order, system_rhs)
end

# module isotropic ends here
end
# module isotropic ends here
