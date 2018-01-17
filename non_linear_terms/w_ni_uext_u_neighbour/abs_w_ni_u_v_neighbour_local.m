function [ res ] = abs_w_ni_u_v_neighbour_local( llcoord, params, paramsP,...
    grid, tria_index, local_vertex_index )
%ABS_W_NI_U_V_NEIGHBOUR_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);

gids = ldg_global_dof_index(params, grid);

N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];

N_neighbour = -N;

tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);

if tria_index_neighbour > 0
    
    local_vertex_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
    
    lcoord_self = llocal2local(grid,local_vertex_index,llcoord);
    
    velocity_basis_self = ldg_evaluate_basis(lcoord_self,params);
    
    lcoord_neighbour = llocal2local(grid,local_vertex_index_neighbour,llcoord);
    
    velocity_basis_neighbour = ldg_evaluate_basis(lcoord_neighbour,params);
    
    w = velocity_basis_neighbour'*params.dofs(gids(tria_index_neighbour,:));
    
    w_dot_N = N_neighbour * w;
    
    w_dot_N = abs(w_dot_N);
    
    res = velocity_basis_self * velocity_basis_neighbour';
    
    res = w_dot_N * res;
    
end

%check 1. return zero for dirichlet boundary

end