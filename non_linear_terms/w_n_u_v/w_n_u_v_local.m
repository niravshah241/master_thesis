function [ res ] = w_n_u_v_local...
    (llcoord, params, paramsP, grid, tria_index, local_vertex_index )
%W_N_U_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);

lcoord = llocal2local(grid,local_vertex_index,llcoord);

N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];

velocity_basis = ldg_evaluate_basis(lcoord,params);

gids = ldg_global_dof_index(params,grid);

dofs = velocity_basis' * params.dofs(gids(tria_index,:));

w_cdot_N = N * dofs;

res = velocity_basis * velocity_basis';

res = w_cdot_N * res;

end