function [res] = error_velocity_rbasis_local(lcoord, params, params_reduced, ...
    B_velocity, grid, tria_index)

velocity_basis = ldg_evaluate_basis(lcoord,params);
gids = ldg_global_dof_index(params, grid);
velocity = velocity_basis' * params.dofs(gids(tria_index,:));
velocity_reduced_basis = velocity_basis' * B_velocity(gids(tria_index,:),:);
velocity_reduced = velocity_reduced_basis * params_reduced.dofs;
res = norm(velocity - velocity_reduced,2)^2;