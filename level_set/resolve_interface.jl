module resolve_interface

using basis, reinitialize_levelset

export interfaceEdgeIntersection

"""
    interfaceEdgeIntersection(distance::Array{Float64, 1},
        fbasis::Basis{T}) where T
return the points ON THE REFERENCE ELEMENT representing the intersection
of the interface and the boundary of the element.
"""
function interfaceEdgeIntersection(distance::Array{Float64, 1},
    fbasis::Basis{T}) where T

    nbr = neighborNodes(T)
    nnodes = length(fbasis.support_points)
    visited = zeros(Bool, nnodes)
    projected = Float64[]
    for I in 1:nnodes
        for J in nbr[:,I]
            if !visited[I] || !visited[J]
                if distance[I]*distance[J] < 0.0
                    XI = fbasis.support_points[I]
                    XJ = fbasis.support_points[J]
                    direction = XJ - XI
                    X0 = 0.5*(XI + XJ)
                    append!(projected, project(X0, distance, direction, fbasis))
                    visited[I] = true
                    visited[J] = true
                end
            end
        end
    end
    return reshape(projected, 2, :)
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
    projected = zeros(size(reference_projected))
    for I in 1:size(projected)[2]
        p = reference_projected[:,I]
        projected[:,I] = interpolate(nodes, p, fbasis)
    end
    return projected
end


# module resolve_interface ends here
end
# module resolve_interface ends here
