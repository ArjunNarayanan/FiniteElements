using WriteVTK, FiniteElements, JLD2


"""
    checkElementSize(slab_length::Float64, slab_width::Float64,
        domain_x::AbstractArray, domain_y::AbstractArray)
Checks if `slab_length` and `slab_width` are contained within `domain_x` and
`domain_y`.
"""
function checkElementSize(slab_length::Float64, slab_width::Float64,
    domain_x::AbstractArray, domain_y::AbstractArray)

    @assert slab_length in domain_x "Slab length must be integer multiple of element size"
    @assert slab_width in domain_y "Slab width must be an integer multiple of element size"
end


"""
    getNodeIDs(row::Int64, col::Int64, nx::Int64)
returns the connectivity of element at position `row, col` in the grid.
"""
function getNodeIDs(row::Int64, col::Int64, nx::Int64)
    bottom_left = (row-1)*nx + col
    bottom_right = bottom_left + 1
    top_left = bottom_left + nx
    top_right = top_left + 1
    return [bottom_left, bottom_right, top_right, top_left]
end

"""
    getConnectivity(nx::Int64, ny::Int64)
returns the connectivity of the 2D elements in the entire grid.
"""
function getConnectivity(nx::Int64, ny::Int64)
    nel_x = nx - 1
    nel_y = ny - 1
    nel = nel_x*nel_y
    count = 1
    conn = zeros(Int64, 4, nel)
    for j in 1:nel_y
        for i in 1:nel_x
            conn[:,count] = getNodeIDs(j, i, nx)
            count += 1
        end
    end
    return conn
end


"""
    getNodes(domain_x::AbstractArray, domain_y::AbstractArray)
returns the nodal coordinates of a rectangular grid defined by the cartesian
product of `domain_x` and `domain_y`.
"""
function getNodes(domain_x::AbstractArray, domain_y::AbstractArray)
    nx = length(domain_x)
    ny = length(domain_y)
    nnodes = nx*ny
    nodes = zeros(2,nnodes)
    count = 1
    for j in 1:ny
        for i in 1:nx
            nodes[:,count] = [domain_x[i], domain_y[j]]
            count += 1
        end
    end
    return nodes
end

"""
    getNodeGroups(nx::Int64, ny::Int64)
returns `bottom, right, top, left` which are arrays with the node IDs on the
bottom, right, top, and left edges.
"""
function getNodeGroups(nx::Int64, ny::Int64)
    nnodes = nx*ny
    bottom = 1:nx
    right = nx:nx:nnodes
    left = 1:nx:(nnodes-nx+1)
    top = (nnodes-nx+1):nnodes
    return bottom, right, top, left
end

"""
    getDomains(slab_length::Float64, slab_width::Float64,
        element_size::Float64)
returns the 1D domains in x and y directions used to construct the structured
mesh.
"""
function getDomains(slab_length::Float64, slab_width::Float64,
    element_size::Float64)

    domain_x = 0.0:element_size:slab_length
    domain_y = 0.0:element_size:slab_width
    return domain_x, domain_y
end

"""
    getElementGroups(nx::Int64, ny::Int64)
return an array with numberinf of all quad elements.
"""
function getElementGroups(nx::Int64, ny::Int64)
    nel_x = nx - 1
    nel_y = ny - 1
    return collect(1:nel_x*nel_y)
end
"""
    generateStructuredMesh(slab_length::Float64, slab_width::Float64,
        element_size::Float64)
returns a `Mesh` object with the required data for the finite element
calculation.
"""
function generateStructuredMesh(slab_length::Float64, slab_width::Float64,
    element_size::Float64)

    domain_x, domain_y = getDomains(slab_length, slab_width, element_size)
    nx = length(domain_x)
    ny = length(domain_y)
    checkElementSize(slab_length, slab_width, domain_x, domain_y)


    data = Dict{Symbol, Any}()

    nodes = getNodes(domain_x::AbstractArray, domain_y::AbstractArray)
    data[:nodes] = nodes

    data[:elements] = Dict{DataType, Array}()

    connectivity = getConnectivity(nx, ny)
    data[:elements][Quadrilateral{4}] = connectivity

    bottom, right, top, left = getNodeGroups(nx, ny)
    data[:node_groups] = Dict{String, Array}()
    data[:node_groups]["bottom"] = bottom
    data[:node_groups]["right"] = right
    data[:node_groups]["top"] = top
    data[:node_groups]["left"] = left

    body_elements = getElementGroups(nx, ny)
    data[:element_groups] = Dict()
    data[:element_groups]["body"] = Dict()
    data[:element_groups]["body"][Quadrilateral{4}] = body_elements

    mesh = Mesh{2}(data)

    return mesh
end


function saveVTUmesh(mesh, filename)
    output = Output(filename, [Quadrilateral{4}], mesh)
    outfiles = vtk_save(output.vtkfile)
end


const slab_length = 1.0
const slab_width = 1.0
const element_size = 1e-1

mesh = generateStructuredMesh(slab_length, slab_width,
    element_size)

nel_x = round(Int, slab_length/element_size)
outfolder = "mesh/Quad/"
try
    mkpath(outfolder)
catch
end


filename = outfolder*"Quad_"*string(nel_x)

# Save a VTU mesh if required
saveVTUmesh(mesh, filename)

# Save a JLD mesh
@save filename*".jld2" mesh
