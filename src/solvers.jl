module solvers

using assembly, IterativeSolvers


function solveDirect(system::GlobalSystem)
	system.D[:] = K\Array(system.F)
end






# module solvers ends here
end
# module solvers ends here