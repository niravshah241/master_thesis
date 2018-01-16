function [ res ] = w_ni_uext_v_local(llcoord, params, paramsP, grid,...
    tria_index, local_vertex_index)
%W_NI_UEXT_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
gids = ldg_global_dof_index(params,grid);
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);

if tria_index_neighbour > 0 
    local_vertex_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
    lcoord_internal = llocal2local(grid,local_vertex_index,llcoord);
    lcoord_external = llocal2local(grid,local_vertex_index_neighbour,llcoord);
    velocity_basis_internal = ldg_evaluate_basis(lcoord_internal,params);
    velocity_basis_external = ldg_evaluate_basis(lcoord_external,params);
    w = velocity_basis_internal' * params.dofs(gids(tria_index,:));
    N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
    w_dot_N = N * w;
    res = velocity_basis_internal * velocity_basis_external';
    res = w_dot_N * res;
end

end
%check 1. % for dirichlet boundary this function returns zeros only 
%2. check formula