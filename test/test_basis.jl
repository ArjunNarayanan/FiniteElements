using Test
using FiniteElements

basis_types = [Line{2}, Line{3}, Triangle{3}, Triangle{6},
				Quadrilateral{4}, Quadrilateral{9}]

tol = 1e-10

for t in basis_types
	basis = Basis(t)
	for i in eachindex(basis.functions)
		f = basis.functions[i]
		for j in eachindex(basis.support_points)
			p = basis.support_points[j]
			if i == j
				@test abs(f(p) - 1.0) < tol
			else
				@test abs(f(p)) < tol
			end
		end
	end
end










true