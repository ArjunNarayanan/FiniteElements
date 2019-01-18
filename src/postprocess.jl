module postprocess
using geometry, WriteVTK

export Output, writePointVectors, writeCellScalars


geom_to_vtk = Dict(Vertex => VTKCellTypes.VTK_VERTEX,
				   Line{2} => VTKCellTypes.VTK_LINE,
				   Line{3} => VTKCellTypes.VTK_QUADRATIC_EDGE,
				   Triangle{3} => VTKCellTypes.VTK_TRIANGLE,
				   Triangle{6} => VTKCellTypes.VTK_QUADRATIC_TRIANGLE,
				   Quadrilateral{4} => VTKCellTypes.VTK_QUAD,
				   Quadrilateral{9} => VTKCellTypes.VTK_QUADRATIC_QUAD)

"""
	Output
Struct to store necessary information to generate output VTK files.
# Attributes
- `spacedim::Int64` - the spatial dimension of the mesh. 
- `vtkfile::WriteVTK.DatasetFile`
- `elTypes::Array{DataType, 1}` - the element types to be used in postprocessing.
# Constructor
	Output(filename::String, elTypes::Array{DataType, 1}, mesh::Mesh{spacedim}) where spacedim
"""
struct Output
	spacedim::Int64
	vtkfile::WriteVTK.DatasetFile
	elTypes::Array{DataType, 1}
	function Output(filename::String, elTypes::Array{DataType, 1}, 
					mesh::Mesh{spacedim}) where spacedim
		points = mesh.data[:nodes]
		cells = Array{MeshCell{Array{Int64, 1}},1}()
		for elType in elTypes
			elements = mesh.data[:elements][elType]
			for e in 1:size(elements)[2]
				cell = MeshCell(geom_to_vtk[elType], elements[:,e])
				push!(cells, cell)
			end
		end
		vtkfile = vtk_grid(filename, points, cells)
		new(spacedim, vtkfile, elTypes)
	end
end


"""
	writePointVectors(output::Output, 
					  vector::Array{Float64, 2},
					  var_name::String)
Add `vector` as point data into `output.vtkfile` with variable name `var_name`.
`vector` must be of size `(dim,num_nodes)` where `dim` is either 2 or 3 and
`num_nodes` is the number of nodes in the mesh used to initialize `output`.
"""
function writePointVectors(output::Output, 
						   vector::Array{Float64, 2},
						   var_name::String)
	if size(vector)[1] == 2
		v_ext = zeros(3, size(vector)[2])
		v_ext[1:2,:] = vector
		vtk_point_data(output.vtkfile, v_ext, var_name)
	elseif size(vector)[1] == 3
		vtk_point_data(output.vtkfile, vector, var_name)
	else
		error("size(vector)[1] should be 2 or 3")
	end
end

"""
	writeCellScalars(output::Output,
					 scalar::Array{Float64, 1},
					 var_name::String)
Add `scalar` as cell data into `output.vtkfile` with variable name `var_name`.
Note that the ordering of `scalar` as cell data is assumed to be consistent with
the ordering used to initialize `output`.
"""
function writeCellScalars(output::Output,
						  scalar::Array{Float64, 1},
						  var_name::String)
	vtk_cell_data(output.vtkfile, scalar, var_name)
end





# module ends here
end
# module ends here