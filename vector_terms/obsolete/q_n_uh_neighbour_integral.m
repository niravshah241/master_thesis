function [ res] = q_n_uh_neighbour_integral( grid, params, tria_index,...
    local_vertex_index, qdeg)
%Q_N_UH_NEIGHBOUR_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
face_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);

f = @(llcoord) q_n_uh_neighbour_local( llcoord, grid, params, tria_index,...
    local_vertex_index);
res = intervalquadrature(qdeg,f)*grid.EL(tria_index_neighbour,...
    face_index_neighbour);

end