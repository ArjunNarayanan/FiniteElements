import gmsh

model = gmsh.model
geo = model.geo

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 1)

model.add("Test")

# Define the side length of the square
a = 1.0
# Define the element size to be used
lc1 = 0.1

# Define the four corners of the square
p1 = geo.addPoint(0.0, 0.0, 0.0, lc1)
p2 = geo.addPoint(a, 0.0, 0.0, lc1)
p3 = geo.addPoint(a, a, 0.0, lc1)
p4 = geo.addPoint(0.0, a, 0.0, lc1)

# Define the line segments of the boundary of the square
l1 = geo.addLine(p1, p2)
l2 = geo.addLine(p2, p3)
l3 = geo.addLine(p3, p4)
l4 = geo.addLine(p4, p1)
# Define the diagonal
d = geo.addLine(p1, p3)

# Define two surfaces
cl1 = geo.addCurveLoop([l1, l2, -d])
cl2 = geo.addCurveLoop([d, l3, l4])
# Define the surface
surf1 = geo.addPlaneSurface([cl1])
surf2 = geo.addPlaneSurface([cl2])



# make a physical group for the surface
pg_surf1 = model.addPhysicalGroup(2, [surf1])
pg_surf2 = model.addPhysicalGroup(2, [surf2])
pg_line = model.addPhysicalGroup(1, [d])
model.setPhysicalName(2, pg_surf1, "surface1")
model.setPhysicalName(2, pg_surf2, "surface2")
model.setPhysicalName(1, pg_line, "line")

geo.synchronize()

model.mesh.setRecombine(2, surf1)
gmsh.option.setNumber("Mesh.Algorithm", 8)

model.mesh.generate(2)

gmsh.write("Test.msh")
gmsh.write("Test.vtk")

gmsh.finalize()


















