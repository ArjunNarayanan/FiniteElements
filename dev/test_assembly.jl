using FiniteElements, StaticArrays

# Define a 3-node linear triangle
p1 = Point(1,1)
p2 = Point(4,2)
p3 = Point(2,4)
t1 = Triangle{1}((1,2,3), (p1,p2,p3))

# Get the master element for this element with second order
# quadrature. Precompute the values and gradients of the basis
master = Master(Triangle{1}, 2, :values, :gradients)

# Compute the mapping for this element
mapping = Map(master, Triangle{1,2}, :coordinates, :derivatives)
reinit(mapping, master, t1)

for q in eachindex(master.quadrature)
	pq, wq = master.quadrature[q]

end