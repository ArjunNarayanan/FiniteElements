using PyCall, geometry
@pyimport meshio


spacedim = 2

elementTypes = Dict("vertex" => Vertex,
					"line" => Line{2},
					"line3" => Line{3},
					"triangle" => Triangle{3},
					"triangle6" => Triangle{6},
					"quad" => Quadrilateral{4},
					"quad9" => Quadrilateral{9})

mesh_data = meshio.read("Test.msh")

mesh = Dict()

mesh[:points] = mesh_data[:points][:,1:spacedim]

for key in keys(mesh_data[:cells])
	eType = elementTypes[key]
	mesh[eType] = mesh_data[:cells][key]
end

tagToKey = Dict()

for key in keys(mesh_data[:field_data])
	tag = mesh_data[:field_data][key][1]
	tagToKey[tag] = key
	mesh[key] = Array{Tuple{UnionAll, Int64}, 1}()
end

for key in keys(mesh_data[:cell_data])
	eType = elementTypes[key]
	tags = mesh_data[:cell_data][key]["gmsh:physical"]
	for i in eachindex(tags)
		tag = tags[i]
		group_name = tagToKey[tag]
		push!(mesh[group_name], (eType, i))
	end
end



