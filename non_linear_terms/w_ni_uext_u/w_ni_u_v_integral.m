function [ res ] = w_ni_u_v_integral( params, paramsP, grid,...
    tria_index, local_vertex_index, qdeg )
%W_NI_U_V_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(llcoord) w_ni_u_v_local( llcoord, params, paramsP, grid,...
    tria_index, local_vertex_index );
res = 1/2*intervalquadrature(qdeg,f)*grid.EL(tria_index,local_vertex_index);


end

