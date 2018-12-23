module elements

using StaticArrays

using geometry, quadrature, master, maps


struct Element
	data::Dict{Symbol, Array}
	r::Ref{Map}
	function Element(mapping::Map)
		r = Ref(mapping)
		for arg in mapping.args
			
		end
	end
end



struct ElementLibrary
	data::Dict{UnionAll, Array}
	function ElementLibrary()
		data = Dict{UnionAll, Array}()
		new(data)
	end
end




# module elements ends here
end
# module elements ends here