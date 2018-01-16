function [ res ] = w_ni_u_v_local( llcoord, params, paramsP, grid,...
    tria_index, local_vertex_index )
%W_NI_U_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);

gids = ldg_global_dof_index(params,grid);

lcoord = llocal2local(grid,local_vertex_index,llcoord);

velocity_basis = ldg_evaluate_basis(lcoord,params);

N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];

w = velocity_basis' * params.dofs(gids(tria_index,:));

w_dot_N = N * w;

res = velocity_basis * velocity_basis';

res = w_dot_N * res;

end