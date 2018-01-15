function [ res ] = c_h_u_v_term_1_local( lcoord, params, paramsP, grid, tria_index)
%C_H_U_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);
res = zeros(params.ndofs_per_element);
velocity_basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);
velocity_basis = ldg_evaluate_basis(lcoord,params);
w = velocity_basis'*params.dofs(gids(tria_index,:));
velocity_basis_vector = sum(velocity_basis,2);

for a=1:1:params.dimrange
    for i=a:params.dimrange:params.ndofs_per_element
        for j=a:params.dimrange:params.ndofs_per_element
            del_v = sum(velocity_basis_derivative{i},2);
            res(i,j) = res(i,j) + w'*del_v*velocity_basis_vector(j);
        end
    end
end
end