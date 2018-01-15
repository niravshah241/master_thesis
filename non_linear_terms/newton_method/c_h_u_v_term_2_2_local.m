function [ res ] = c_h_u_v_term_2_2_local...
    ( llcoord, params, paramsP, grid, tria_index, local_vertex_index)
%C_H_U_V_TERM_2_1_LOCAL Summary of this function goes here
%   Detailed explanation goes here
res = zeros(params.ndofs_per_element);
% checkpoint 1: currently assumes u_ext for dirichlet boundyr zero
gids = ldg_global_dof_index(params,grid);

lcoord = llocal2local(grid,local_vertex_index,llcoord);
velocity_basis = ldg_evaluate_basis(lcoord,params);

w = velocity_basis'*params.dofs(gids(tria_index,:));%checkpoint 2
N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
w_dot_n = N*w;

res = velocity_basis * velocity_basis';
res = w_dot_n * res;

end