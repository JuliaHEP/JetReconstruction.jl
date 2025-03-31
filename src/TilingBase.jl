#using JetReconstruction

abstract type TilingBase end

tile_index(p::PseudoJet)::Int64 = begin 
    return 0
end 

n_tiles()::Int64 = begin 
   return 0
end 

n_good_tiles()::Int64 = begin 
    return n_tiles()
end 

tile_is_good(itile::Int64)::Bool = begin 
    return true
end 

all_tiles_good()::Bool = begin
    return n_good_tiles() == n_tiles()
end 

all_tiles_equal_area()::Bool =begin
    return true
end 

tile_area(itile::Int64)::Float64 = begin
    return mean_tile_area()
end

mean_tile_area()::Float64 = begin 
    return 0
end 

description()::String = begin 
    return 0
end  

is_initialised():: String = begin 
    return 0
end  

is_initialized()::Bool = begin
    return is_initialised()
end 



