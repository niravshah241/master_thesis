function [ res ] = w_ni_uext_v_neighbour_local( llcoord, params, paramsP,...
    grid, tria_index, local_vertex_index)
%W_NI_UEST_V_NEIGHBOUR_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);

gids = ldg_global_dof_index(params, grid);

tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);

if tria_index_neighbour > 0

local_vertex_index_neighbour =...
    find(grid.NBI(tria_index_neighbour,:)==tria_index);
    
tria_index_neighbour_ext =...
    grid.NBI(tria_index_neighbour,local_vertex_index_neighbour);

local_vertex_index_neighbour_ext =...
    find(grid.NBI(tria_index_neighbour_ext,:)==tria_index_neighbour);

if tria_index_neighbour_ext > 0

N = [grid.NX(tria_index,local_vertex_index) ...
    grid.NY(tria_index,local_vertex_index)];

N_neighbour = -N;

N_neighbour_ext = -N_neighbour;

lcoord = llocal2local(grid,local_vertex_index,llcoord);
lcoord_neighbour = llocal2local(grid,local_vertex_index_neighbour,llcoord);
lcoord_neighbour_ext = llocal2local(grid,local_vertex_index_neighbour_ext,llcoord);

velocity_basis_neighbour_ext = ldg_evaluate_basis(lcoord_neighbour_ext,params);

dofs_ext = velocity_basis_neighbour_ext'...
    *params.dofs(gids(tria_index_neighbour_ext,:));

w_dot_N = N_neighbour_ext * dofs_ext;

velocity_basis_self = ldg_evaluate_basis(lcoord,params);

res = velocity_basis_self * velocity_basis_neighbour_ext';

res = w_dot_N * res;

end
end
end