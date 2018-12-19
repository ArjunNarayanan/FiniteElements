using Test
using geometry, elements

basis_types = [Line{1}, Line{2}, Triangle{1}, Triangle{2},
				Quadrilateral{1}, Quadrilateral{2}]

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