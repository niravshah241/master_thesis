function [ res ] = non_linear_term_neighbour_assembly( params,paramsP,grid,qdeg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res_internal = zeros(params.ndofs);
res_external = zeros(params.ndofs);
res_internal_abs = zeros(params.ndofs);
res_external_abs = zeros(params.ndofs);

gids = ldg_global_dof_index(params,grid);

[ tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet,...
    b_index_dirichlet, local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

[ tria_index_internal,local_edge_index_internal,a_index_internal,...
    b_index_internal, local_vertex_index_internal] = tria_edge_index_internal( grid );

tria_index = [tria_index_internal,tria_index_dirichlet];
local_vertex_index = [local_vertex_index_internal,local_vertex_index_dirichlet];

for i = 1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_internal(gids(tria_index(i),:),gids(tria_index_neighbour,:)) = ...
            res_internal(gids(tria_index(i),:),gids(tria_index_neighbour,:)) + ...
            w_ni_u_v_neighbour_integral( params, paramsP, grid, tria_index(i),...
            local_vertex_index(i), qdeg );
        res_internal_abs(gids(tria_index(i),:),gids(tria_index_neighbour,:)) = ...
            res_internal_abs(gids(tria_index(i),:),gids(tria_index_neighbour,:)) + ...
            abs_w_ni_u_v_neighbour_integral(params, paramsP, grid,...
            tria_index(i), local_vertex_index(i), qdeg);
        local_vertex_index_neighbour = ...
            find(grid.NBI(tria_index_neighbour,:)==tria_index(i));
        tria_index_neighbour_ext = ...
            grid.NBI(tria_index_neighbour,local_vertex_index_neighbour);
        local_vertex_index_neighbour_ext = ...
            find(grid.NBI(tria_index_neighbour_ext,:)== local_vertex_index_neighbour);
        if tria_index_neighbour_ext > 0
            res_external(gids(tria_index_neighbour,:),gids(tria_index_neighbour_ext,:)) = ...
                res_external(gids(tria_index_neighbour,:),gids(tria_index_neighbour_ext,:)) + ...
                w_ni_uext_v_neighbour_integral( params, paramsP, grid,...
                tria_index(i), local_vertex_index(i), qdeg);
            res_external_abs(gids(tria_index_neighbour,:),gids(tria_index_neighbour_ext,:)) = ...
                res_external(gids(tria_index_neighbour,:),gids(tria_index_neighbour_ext,:)) + ...
                abs_w_ni_uext_v_neighbour_integral( params, paramsP, grid,...
                tria_index(i), local_vertex_index(i), qdeg);
        end
    end
end

res.res_internal = zeros(params.ndofs);
res.res_external = zeros(params.ndofs);
res.res_internal_abs = zeros(params.ndofs);
res.res_external_abs = zeros(params.ndofs);
res.res = zeros(params.ndofs);

res.res_internal = sparse(res_internal);
res.res_external = sparse(res_external);
res.res_internal_abs = sparse(res_internal_abs);
res.res_external_abs = sparse(res_external_abs);
res.res = sparse(res.res_internal + res.res_external +...
    res.res_internal_abs - res.res_external_abs);

if params.show_sparsity == true
    figure()
    spy(full(res.res))
    title('spy of non linear term c(u;u,\phi (neighbour part))')
end

% check assembly indexes

end