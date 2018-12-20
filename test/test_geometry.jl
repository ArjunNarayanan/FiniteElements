using Revise
using PyCall
using FiniteElements, Test
@pyimport meshio

# Test the type parametrization
p1 = Point(0.5)
@test typeof(p1) <: Point{1}

p2 = Point(0.5, 0.6)
@test typeof(p2) <: Point{2}

p3 = Point(1.0, 2.5, 3.5)
@test typeof(p3) <: Point{3}

p4 = Point((1.0, 2.5, 3.5))
@test typeof(p3) <: Point{3}

# Test multiplication by scalar
p4 = 2.0*p3
@test p4.x == [2.0, 5.0, 7.0]
p4 = p3*2.0
@test p4.x == [2.0, 5.0, 7.0]

# Test addition
p5 = p4 + p3
@test p5.x == [3.0, 7.5, 10.5]

# Test subtraction
p6 = p4 - p3
@test p6.x == [1.0, 2.5, 3.5]

# Create a Vertex in 1D, 2D
v1 = Vertex((1,), (p1,))
@test typeof(v1) <: Vertex{0, 1}
v1 = Vertex((1,), (p2,))
@test typeof(v1) <: Vertex{0, 2}

# Create line elements in 1D and 2D
p1 = Point(0.5)
p2 = Point(2.0)
l1 = Line((1,2), (p1, p2))
@test typeof(l1) <: Line{1, 1}
@test length(l1.nodes) == 2

p3 = Point(3.0)
l2 = Line((1,2,3), (p1,p2,p3))
@test typeof(l2) <: Line{2, 1}
@test length(l2.nodes) == 3

p1 = Point(0.0, 0.0)
p2 = Point(1.0, 0.0)
l1 = Line((1,2), (p1, p2))
@test typeof(l1) <: Line{1, 2}
@test length(l1.nodes) == 2

p3 = Point(0.5, 0.0)
l2 = Line((1,2,3), (p1, p2, p3))
@test typeof(l2) <: Line{2, 2}
@test length(l2.nodes) == 3

# Create triangle elements in 2D
p1 = Point(0.0, 0.0)
p2 = Point(1.0, 0.0)
p3 = Point(0.0, 1.0)
t1 = Triangle((1,2,3), (p1,p2,p3))
@test typeof(t1) <: Triangle{1, 2}
@test length(t1.nodes) == 3


p1 = Point(0.0, 0.0)
p2 = Point(1.0, 0.0)
p3 = Point(0.0, 1.0)
p4 = 0.5*(p1 + p2)
p5 = 0.5*(p2 + p3)
p6 = 0.5*(p1 + p3)
t2 = Triangle((1,2,3,4,5,6),
				(p1,p2,p3,p4,p5,p6))
@test typeof(t2) <: Triangle{2, 2}
@test length(t2.nodes) == 6

# Create Quad elements in 2D
p1 = Point(0.0, 0.0)
p2 = Point(1.0, 0.0)
p3 = Point(1.0, 1.0)
p4 = Point(0.0, 1.0)
q1 = Quadrilateral((1,2,3,4),
			(p1,p2,p3,p4))
@test typeof(q1) <: Quadrilateral{1, 2}
@test length(q1.nodes) == 4

p1 = Point(0.0, 0.0)
p2 = Point(1.0, 0.0)
p3 = Point(1.0, 1.0)
p4 = Point(0.0, 1.0)
p5 = 0.5*(p1 + p2)
p6 = 0.5*(p2 + p3)
p7 = 0.5*(p3 + p4)
p8 = 0.5*(p4 + p1)
p9 = Point(0.5, 0.5)
q2 = Quadrilateral((1,2,3,4,5,6,7,8,9),
			(p1,p2,p3,p4,p5,p6,p7,p8,p9))
@test typeof(q2) <: Quadrilateral{2, 2}
@test length(q2.nodes) == 9

mesh_data = meshio.read("Test.msh")
spacedim = 2
mesh = LoadMesh(mesh_data, spacedim)
nElements = sum([size(mesh_data[:cells][key])[1] for key in keys(mesh_data[:cells])])

p_rand_index = rand(1:size(mesh_data[:points])[1])
prand = Point(tuple([mesh_data[:points][p_rand_index,j] for j in 1:spacedim]...))

@test length(mesh.elements) == nElements
@test prand == mesh.nodes[p_rand_index]

@test haskey(mesh.element_groups, "surface1")
@test haskey(mesh.element_groups, "surface2")
@test haskey(mesh.element_groups, "line")

# @test haskey(mesh.element_types,)




# If all tests pass
true