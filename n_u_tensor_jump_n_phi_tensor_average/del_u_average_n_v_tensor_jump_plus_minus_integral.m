function [ res ] = del_u_average_n_v_tensor_jump_plus_minus_integral...
    (params,grid,tria_index,local_vertex_index,qdeg)

f = @(llcoord) del_u_average_n_v_tensor_jump_plus_minus_local( llcoord,...
    params,grid,tria_index,local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index)/2;

%grid.EL(tria_index,local_vertex_index) = ...
%grid.EL(tria_index_neighbour,local_vertex_index_neighbour);