function [ res ] = w_n_u_v_assembly( params, paramsP, grid, qdeg)
%W_N_U_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs);

[tria_index_neumann,local_edge_index_neumann,a_index_neumann,...
    b_index_neumann, local_vertex_index_neumann] = tria_edge_index_neumann( grid );


gids = ldg_global_dof_index(params,grid);

for i = 1:1:length(tria_index_neumann)
    res(gids(tria_index_neumann(i),:),gids(tria_index_neumann(i),:))...
        = res(gids(tria_index_neumann(i),:),gids(tria_index_neumann(i),:))...
        + w_n_u_v_integral(params, paramsP, grid, tria_index_neumann(i),...
        local_vertex_index_neumann(i), qdeg);
end

if params.show_sparsity == true
    figure()
    spy(res)
    title('(w \cdot n) u \cdot \phi')
end

end