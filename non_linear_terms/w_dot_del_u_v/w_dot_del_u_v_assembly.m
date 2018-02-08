function [ res ] = w_dot_del_u_v_assembly(params, paramsP, grid, qdeg )
%W_DOT_DEL_U_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs);
gids = ldg_global_dof_index(params,grid);

for tria_index = 1:1:grid.nelements
    res(gids(tria_index,:),gids(tria_index,:)) =...
        res(gids(tria_index,:),gids(tria_index,:)) + ...
        w_dot_del_u_v_integral( params, paramsP, grid, tria_index, qdeg);
end

if params.show_sparsity == true
    spy(res)
    title('((u_k \cdot \nabla)\phi,\phi)')
    axis tight
    axis equal
end

end

