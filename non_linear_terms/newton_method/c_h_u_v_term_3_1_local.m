function [ res ] = c_h_u_v_term_3_1_local...
    ( llcoord, params, paramsP, grid, tria_index, local_vertex_index)
%C_H_U_V_TERM_2_1_LOCAL Summary of this function goes here
%   Detailed explanation goes here
res = zeros(params.ndofs_per_element);
% checkpoint 1: currently assumes u_ext for dirichlet boundyr zero
gids = ldg_global_dof_index(params,grid);
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
if tria_index_neighbour>0
    local_vertex_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
    lcoord = llocal2local(grid,local_vertex_index,llcoord);
    velocity_basis = ldg_evaluate_basis(lcoord,params);
    
    lcoord_neighbour = llocal2local(grid,local_vertex_index_neighbour,llcoord);
    velocity_basis_neighbour = ldg_evaluate_basis(lcoord_neighbour,params);
    
    w = velocity_basis'*params.dofs(gids(tria_index,:));%checkpoint 2
    N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
    w_dot_n = abs(N*w);
    
    res = velocity_basis * velocity_basis_neighbour';
    res = w_dot_n * res;
    
end