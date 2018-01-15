function [ res ] = w_dot_del_u_v_local( lcoord, params, paramsP, grid, tria_index)
%W_DOT_DEL_U_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);

gids = ldg_global_dof_index(params,grid);

velocity_basis = ldg_evaluate_basis(lcoord,params);

dofs = velocity_basis' * params.dofs(gids(tria_index,:));

velocity_basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);

velocity_basis_scalar = sum(velocity_basis,2);

for a=1:1:params.dimrange
    for i=a:params.dimrange:params.ndofs_per_element
        temp = (sum(velocity_basis_derivative{i},1))*dofs;
        for j=a:params.dimrange:params.ndofs_per_element
            res(i,j) = res(i,j) + temp * velocity_basis_scalar(j);
        end
    end
end

end