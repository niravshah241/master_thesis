function [ res ] = c_h_u_v_term_4_local...
    ( llcoord, params, paramsP, grid, tria_index, local_vertex_index )
%C_H_U_V_TERM_4 Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
gids = ldg_global_dof_index(params,grid);
lcoord = llocal2local(grid,local_vertex_index,llcoord);
velocity_basis = ldg_evaluate_basis(lcoord,params);
N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
w = ldg_evaluate_basis(lcoord,params)'*params.dofs(gids(tria_index,:));
w_dot_n = N*w;
for i=1:1:params.ndofs_per_element
    for j=1:1:params.ndofs_per_element
        res(i,j) = res(i,j) + velocity_basis(i,:)*velocity_basis(j,:)';
    end
end

res = res*w_dot_n;

end