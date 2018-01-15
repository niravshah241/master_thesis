function [ res ] = c_h_u_v_term_4_integral...
    ( params, paramsP, grid, tria_index, local_vertex_index, qdeg)
%C_H_U_V_TERM_4_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(llcoord)  c_h_u_v_term_4_local...
    ( llcoord, params, paramsP, grid, tria_index, local_vertex_index );
res = intervalquadrature(qdeg,f)*grid.EL(tria_index,local_vertex_index);

end

