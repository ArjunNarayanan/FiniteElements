using Revise
using PyCall
using geometry, Test
@pyimport meshio





@test Vertex <: Triangulation{1, 0}

@test Line{2} <: Triangulation{2, 1}

@test Line{3} <: Triangulation{3, 1}

@test Triangle{3} <: Triangulation{3, 2}

@test Triangle{6} <: Triangulation{6, 2}

@test Quadrilateral{4} <: Triangulation{4, 2}

@test Quadrilateral{9} <: Triangulation{9, 2}




mesh_data = meshio.read("Test.msh")
spacedim = 2
mesh = Mesh{2}(mesh_data)

nElements_expected = sum([size(mesh_data[:cells][key])[1] for key in keys(mesh_data[:cells])])
nElements = sum([size(mesh.data[:cells][key])[1] for key in keys(mesh.data[:cells])])

@test nElements == nElements_expected



@test haskey(mesh.data[:groups], "surface1")
@test haskey(mesh.data[:groups], "surface2")
@test haskey(mesh.data[:groups], "line")

@test haskey(mesh.data[:cells], Triangle{3})
@test haskey(mesh.data[:cells], Line{2})
@test haskey(mesh.data[:cells], Quadrilateral{4})




# If all tests pass
true