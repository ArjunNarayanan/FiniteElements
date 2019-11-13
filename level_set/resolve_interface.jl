module resolve_interface

using basis, reinitialize_levelset, LinearAlgebra

export interfaceEdgeIntersection, fitNormal

"""
    interfaceEdgeIntersection(distance::Array{Float64, 1},
        fbasis::Basis{T}) where T
return a point ON THE REFERENCE ELEMENT representing the intersection
of the interface and the boundary of the element.
"""
function interfaceEdgeIntersection(distance::Array{Float64, 1},
    fbasis::Basis)

    nnodes = length(fbasis.support_points)
    neighbors = circshift(1:nnodes, -1)

    projected = zeros(length(fbasis.support_points[1]))
    for node1 in 1:nnodes
        node2 = neighbors[node1]
        if distance[node1]*distance[node2] <= 0.0
            X1 = fbasis.support_points[node1]
            X2 = fbasis.support_points[node2]
            D1 = distance[node1]
            D2 = distance[node2]
            direction = (X2 - X1)/norm(X2 - X1)
            X0 = 0.5*(X1 + X2)
            projected = project(X0, distance, direction, fbasis)
        end
    end
    return projected
end

"""
    interfaceEdgeIntersection(nodes::Array{Float64, 2},
        distance::Array{Float64, 1}, fbasis::Basis)
return the points on the element defined by `nodes` representing the
intersection of the interface and the boundary of the element.
"""
function interfaceEdgeIntersection(nodes::Array{Float64, 2},
    distance::Array{Float64, 1}, fbasis::Basis)

    reference_projected = interfaceEdgeIntersection(distance, fbasis)
    projected = interpolate(nodes, reference_projected, fbasis)
    return projected
end



"""
    fitNormal(nodes::Array{Float64, 2}, distance::Array{Float64, 1},
        x0::Array{Float64, 1})
returns the `normal` to a hyperplane passing through `x0` such that the
distance of `nodes` to the hyperplane is approximately `distance` (in a least-
squares sense).
"""
function fitNormal(nodes::Array{Float64, 2}, distance::Array{Float64, 1},
    x0::Array{Float64, 1})

    A = nodes .- x0
    return (A*A')\(A*distance)
end


# module resolve_interface ends here
end
# module resolve_interface ends here
