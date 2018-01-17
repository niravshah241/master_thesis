function [ res ] = w_ni_uext_v_neighbour_integral( params, paramsP,...
    grid, tria_index, local_vertex_index, qdeg)
%W_NI_UEXT_V_NEIGHBOUR_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here


f = @(llcoord) w_ni_uext_v_neighbour_local( llcoord, params, paramsP,...
    grid, tria_index, local_vertex_index);

res = 1/2 * intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);

% check: 1. intergal is taken on grid.EL(tria_index,local_vertex_index), or
% grid.EL(tria_index_neighbour,local_vertex_index_neighbour) or 
% grid.EL(tria_index_neighbour_ext,local_vertex_index_neighbour_ext)

end

