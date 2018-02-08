function [ res ] = abs_w_ni_u_v_assembly( params, paramsP, grid, qdeg )
%W_NI_U_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res_internal = zeros(params.ndofs);
res_external = zeros(params.ndofs);
gids = ldg_global_dof_index(params,grid);

[ tria_index_dirichlet, local_edge_index_dirichlet, a_index_dirichlet,...
    b_index_dirichlet, local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

[ tria_index_internal, local_edge_index_internal, a_index_internal,...
    b_index_internal, local_vertex_index_internal] = tria_edge_index_internal( grid );


tria_index = [tria_index_internal,tria_index_dirichlet];
local_vertex_index = [local_vertex_index_internal,local_vertex_index_dirichlet];

for i = 1:1:length(tria_index)
    res_internal(gids(tria_index(i),:),gids(tria_index(i),:)) =...
        res_internal(gids(tria_index(i),:),gids(tria_index(i),:)) + ...
        abs_w_ni_u_v_integral( params, paramsP, grid, tria_index(i),...
        local_vertex_index(i), qdeg );
end

for i = 1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_external(gids(tria_index(i),:),gids(tria_index_neighbour,:)) = ...
            res_external(gids(tria_index(i),:),gids(tria_index_neighbour,:)) + ...
            abs_w_ni_uext_v_integral( params, paramsP, grid, tria_index(i),...
            local_vertex_index(i), qdeg);
    end
end

res.res_internal = sparse(res_internal);
res.res_external = sparse(res_external);
res.res = sparse(res_internal - res_external);

if params.show_sparsity == true
    figure()
    spy(full(res_internal))
    title('spy of abs(u_k \cdot n) \phi \cdot \phi')
    figure()
    spy(full(res_internal))
    title('spy of abs(u_k \cdot n) {\phi}^{ext} \cdot \phi')
    figure()
    spy(full(res.res))
    title('spy of abs(u_k \cdot n) \phi \cdot \phi + abs(u_k \cdot n) \phi^{ext} \cdot \phi')
end

end
%check 1. formula and 2. assembly indexes