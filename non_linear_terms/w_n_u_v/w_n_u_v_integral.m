function [ res] = w_n_u_v_integral...
    (params, paramsP, grid, tria_index, local_vertex_index, qdeg)
%W_N_U_V_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(llcoord) w_n_u_v_local(llcoord, params, paramsP,...
    grid, tria_index, local_vertex_index );

res = intervalquadrature(qdeg,f)*grid.EL(tria_index,local_vertex_index);

end

