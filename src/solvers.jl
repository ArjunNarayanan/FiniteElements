module solvers

using assembly

export solveDirect


function solveDirect(system::GlobalSystem)
	system.D[:] = system.K\Array(system.F)
end






# module solvers ends here
end
# module solvers ends here
