using geometry, Test





@test Vertex <: Triangulation{1, 0}

@test Line{2} <: Triangulation{2, 1}

@test Line{3} <: Triangulation{3, 1}

@test Triangle{3} <: Triangulation{3, 2}

@test Triangle{6} <: Triangulation{6, 2}

@test Quadrilateral{4} <: Triangulation{4, 2}

@test Quadrilateral{9} <: Triangulation{9, 2}







# If all tests pass
true