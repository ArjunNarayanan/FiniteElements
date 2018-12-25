using FiniteElements, Test, TensorOperations

const dofs = 2
# Construct an arbitrary Triangle{1} object
p1 = Point(1,1)
p2 = Point(4,2)
p3 = Point(2,4)
T = Triangle{3}((1,2,3), (p1,p2,p3))


element_matrix = elementMatrix(typeof(T), dofs)
element_rhs = elementRHS(typeof(T), dofs)

# Test that the element matrix is of correct size
@test size(element_matrix) == (3,3)
@test size(element_rhs) == (3,)

# Check that each node has the right sized matrix and vector
for m in element_matrix
	@test size(m) == (dofs,dofs)
end
for v in element_rhs
	@test size(v) == (dofs,)
end


# Check reinitialization
for m in element_matrix
	fill!(m, 1.0)
end
for v in element_rhs
	fill!(v, 1.0)
end
reinit(element_matrix)
reinit(element_rhs)

for m in element_matrix
	@test m == zeros(dofs,dofs)
end
for v in element_rhs
	@test v == zeros(dofs)
end



# Now let us try assembly
reinit(element_matrix)
reinit(element_rhs)

master = Master(Triangle{3}, 2, 0, 1)
mapping = Map(master, 2, :coordinates, :derivatives)
reinit(mapping, T)

system_matrix = SystemMatrix(dofs)
system_rhs = SystemRHS(dofs)
KIJ = zeros(2,2)

for q in eachindex(master.quadrature.points)
	pq = master.quadrature.points[q]
	wq = master.quadrature.weights[q]
	for I in 1:length(T.nodes)
		∇ϕI = mapping[:inverse_jacobian][q]'*master[1][I,q]
		for J in 1:length(T.nodes)
			∇ϕJ = mapping[:inverse_jacobian][q]'*master[1][J,q]
			fill!(KIJ, 0.0)
			@tensor begin
				KIJ[i,j] = ∇ϕI[i]*∇ϕJ[j]
			end
			element_matrix[I,J] += KIJ*mapping[:determinant]*wq
		end
		element_rhs[I] += mapping[:inverse_jacobian]'*
							master[1]*mapping[:determinant]*wq
	end
end

updateSystemMatrix(system_matrix, element_matrix)
updateSystemRHS(system_rhs, element_rhs)



true