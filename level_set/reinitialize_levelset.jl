module reinitialize_levelset

using LinearAlgebra, NearestNeighbors
using FiniteElements

export stepProjection, project, hasInterface, projectCentroid,
	distanceToHyperplane, projectCentroidAndComputeNormal, signedDistance,
	nearestNeighbor

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
    basis::Basis{T}; tol = 1e-3, maxiter = 200) where T

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
    end
	return xnext
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
`true` they are of opposite sign or if maximum of `values` is zero.
"""
function hasInterface(values::Array{Float64, 1})
	max_val = maximum(values)
	min_val = minimum(values)
    if max_val*min_val <= 0 && min_val < 0
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
    projectCentroidAndComputeNormal(nodes::Array{Float64, 2},
		values::Array{Float64, 1}, basis::Basis{T}) where T
project the centroid of the reference element of type `T` onto the zero level
set of the interpolation of `values` by `basis`. Compute the normal to the
zero level set at this point. Interpolate the point and the normal into
spatial coordinates. Return the spatial point and normal.
"""
function projectCentroidAndComputeNormal(nodes::Array{Float64, 2},
	values::Array{Float64, 1}, basis::Basis{T}) where T

    reference_centroid = centroid(T)
    projected_centroid = project(reference_centroid, values, basis)
	reference_normal = gradient(values, projected_centroid, basis)
	jacobian = gradient(nodes, projected_centroid, basis)
	spatial_normal = (jacobian')\reference_normal
	spatial_normal = spatial_normal/norm(spatial_normal)

    projected_centroid_spatial = interpolate(nodes, projected_centroid, basis)
    return projected_centroid_spatial, spatial_normal
end

"""
	nearestNeighbor(search_points::Array{Float64, 2},
		query_points::Array{Float64, 2})
wraps the `knn` function from `NearestNeighbors`. Returns the index of the
nearest point in `search_points` to each point in `query_points`.
"""
function nearestNeighbor(search_points::Array{Float64, 2},
	query_points::Array{Float64, 2})

	tree = KDTree(search_points)
	nearest_idx, nearest_dist = knn(tree, query_points, 1)
	nearest_idx = [a[1] for a in nearest_idx]
	return nearest_idx
end

"""
	signedDistance(nodes::Array{Float64, 2}, point_cloud::Array{Float64, 2},
		normals::Array{Float64, 2}, nearest_idx::Array{Int64, 1})
compute the signed distance from each `nodes[:,i]` to a hyperplane passing
through `point_cloud[:,nearest_idx[i]]` with normal `normals[:,nearest_idx[i]]`,
where `nearest_idx[i]` is the index of the nearest point in `point_cloud` to
`nodes[:,i]`.
return the vector of signed distances.
"""
function signedDistance(nodes::Array{Float64, 2}, point_cloud::Array{Float64, 2},
	normals::Array{Float64, 2}, nearest_idx::Array{Int64, 1})

	num_nodes = size(nodes)[2]
	distance = zeros(num_nodes)
	for i in 1:num_nodes
		nearest_point = point_cloud[:,nearest_idx[i]]
		nearest_point_normal = normals[:,nearest_idx[i]]
		distance[i] = distanceToHyperplane(nodes[:,i], nearest_point,
			nearest_point_normal)
	end
	return distance
end

# reinitialize_levelset module ends here
end
# reinitialize_levelset module ends here
