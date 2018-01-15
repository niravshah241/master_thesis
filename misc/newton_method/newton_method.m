function [ res ] = newton_method(params, paramsP, grid, iter_max,...
    tol_solution, tol_update)
%NEWTON_FUNC Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);

for i=1:1:iter_max
    old_solution = params.dofs;
    for k = 1:1:grid.nelements
        func = newton_func(params,gids,k);
        gradient = newton_gradient(params,gids,k);
        params.dofs(gids(k,:)) = params.dofs(gids(k,:)) -...
            inv(gradient) * func;
    end
    if norm(old_solution-params.dofs,2)<tol_update ||...
            norm(params.dofs,2)<tol_solution
        break
    end
end
end