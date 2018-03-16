function [ res ] = error_h0_norm_local( lcoord, params, grid, tria_index)
%ERROR_H0_NORM_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = 0;
JIT = [grid.JIT(tria_index,:,1)',grid.JIT(tria_index,:,2)']; % Jacobian Inverse Transpose
gids = ldg_global_dof_index(params,grid);
basis = ldg_evaluate_basis(lcoord,params);
basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);
dofs = params.dofs(gids(tria_index,:));
solution_derivative = zeros(params.dimrange,2); % 2 is for triagrid
glob = local2global(grid,tria_index,lcoord,params);

for i = 1:1:params.ndofs_per_element
    solution_derivative = solution_derivative + (basis_derivative{i}*JIT'*dofs(i));
end

error = solution_derivative - params.dof_derivative_analytical(glob);


for j = 1:1:params.dimrange
    res = res + norm(error(j,:),2)^2;
end

end