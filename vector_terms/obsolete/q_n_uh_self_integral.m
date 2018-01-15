function [ res] = q_n_uh_self_integral( grid, params, tria_index,...
    local_vertex_index, qdeg)
%Q_N_UH_SELF_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(llcoord) q_n_uh_self_local( llcoord, grid, params, tria_index,...
    local_vertex_index);
res = intervalquadrature(qdeg,f)*grid.EL(tria_index,local_vertex_index);

end

