module reinitialize_levelset

using basis, LinearAlgebra

export stepProjection, project, hasInterface, projectCentroid,
	distanceToHyperplane

"""
    stepProjection(xi::Array{Float64, 1}, values::Array{Float64, 1},
        basis::Basis)
perform one step of Newton interation in projecting `xi` onto the
zero level set of the interpolation of `values` by `basis`.
"""
function stepProjection(xi::Array{Float64, 1}, values::Array{Float64, 1},
    basis::Basis)

    phi = interpolate(values, xi, basis)
    gradphi = gradient(values, xi, basis)
    return xi - phi/norm(gradphi)*gradphi
end

"""
    project(x0::Array{Float64, 1}, values::Array{Float64, 1},
        basis::Basis{T}; tol = 1e-3, maxiter = 100) where T
orthogonally project the point `x0` onto the zero level set of the
interpolation of `values` by `basis` using Newton iterations. The
point `x0` is ASSUMED TO BE IN THE REFERENCE ELEMENT COORDINATES.
"""
function project(x0::Array{Float64, 1}, values::Array{Float64, 1},
    basis::Basis{T}; tol = 1e-3, maxiter = 20) where T

    elDia = diameter(T)
    xprev = copy(x0)
    xnext = stepProjection(xprev, values, basis)
    err = norm(xnext - xprev)/elDia
    count = 1
    while err > tol && count <= maxiter
        xprev = xnext
        xnext = stepProjection(xprev, values, basis)
        err = norm(xnext - xprev)/elDia
        count += 1
    end
    if count > maxiter
        error("Newton iteration failed to converge: $count iterations,
            $err relative error")
    else
        return xnext
    end
end

"""
    stepProjection(xi::Array{Float64, 1}, values::Array{Float64, 1},
        basis::Basis)
perform one step of Newton interation in projecting `xi` onto the
zero level set of the interpolation of `values` by `basis`.
"""
function stepProjection(xi::Array{Float64, 1}, values::Array{Float64, 1},
    direction::Array{Float64, 1}, basis::Basis)

    phi = interpolate(values, xi, basis)
    gradphi = gradient(values, xi, basis)
    return xi - phi/(gradphi'*direction)*direction
end

"""
    project(x0::Array{Float64, 1}, values::Array{Float64, 1},
        direction::Array{Float64, 1}, basis::Basis{T};
        tol = 1e-10, maxiter = 20) where T
project the point `x0` onto the zero level set of the interpolation of `values`
along `direction` by `basis` using Newton iterations.
`direction` is expected to be normalized.
The point `x0` is ASSUMED TO BE IN THE REFERENCE ELEMENT COORDINATES.
"""
function project(x0::Array{Float64, 1}, values::Array{Float64, 1},
    direction::Array{Float64, 1}, basis::Basis{T}; tol = 1e-10,
    maxiter = 20) where T

    elDia = diameter(T)
    xprev = copy(x0)
    xnext = stepProjection(xprev, values, direction, basis)
    err = norm(xnext - xprev)/elDia
    count = 1
    while err > tol && count <= maxiter
        xprev = xnext
        xnext = stepProjection(xprev, values, direction, basis)
        err = norm(xnext - xprev)/elDia
        count += 1
    end
    if count > maxiter
        error("Newton iteration failed to converge: $count iterations,
            $err relative error")
    else
        return xnext
    end
end


"""
    hasInterface(values::Array{Float64, 1})
returns `false` if the maximum and minimum of `values` have the same sign,
`true` if either are zero or of opposite signs.
"""
function hasInterface(values::Array{Float64, 1})
    p = minimum(values)*maximum(values)
    if p <= 0
        return true
    else
        return false
    end
end

"""
	distanceToHyperplane(node::Array{Float64, 1}, x0::Array{Float64, 1},
		normal::Array{Float64, 1})
compute the (signed) distance of `node` to the hyperplane passing through
`x0` with normal `normal`.
"""
function distanceToHyperplane(node::Array{Float64, 1}, x0::Array{Float64, 1},
	normal::Array{Float64, 1})

	return (node - x0)'*normal
end

"""
	distanceToHyperplane(nodes::Array{Float64, 2}, x0::Array{Float64, 1},
		normal::Array{Float64, 1})
calls `distanceToHyperplane` for each vertex in `nodes`.
"""
function distanceToHyperplane(nodes::Array{Float64, 2}, x0::Array{Float64, 1},
	normal::Array{Float64, 1})

	nnodes = size(nodes)[2]
	distance = [distanceToHyperplane(nodes[:,i], x0, normal) for i in 1:nnodes]
	return distance
end

"""
    projectCentroid(nodes::Array{Float64, 2}, values::Array{Float64, 1},
        basis::Basis{T}) where T
project the centroid of the reference element of type `T` onto the zero level
set of the interpolation of `values` by `basis`.
"""
function projectCentroid(nodes::Array{Float64, 2}, values::Array{Float64, 1},
    basis::Basis{T}) where T

    ξ0 = centroid(T)
    ξp = project(ξ0, values, basis)
    xp = interpolate(nodes, ξp, basis)
    return xp
end


# reinitialize_levelset module ends here
end
# reinitialize_levelset module ends here
