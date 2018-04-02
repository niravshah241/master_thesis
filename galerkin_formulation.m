function [params_reduced, paramsP_reduced, stifness_matrix_reduced, rhs_reduced] = ...
    galerkin_formulation(stifness_matrix, rhs, params, paramsP, grid, ...
    red_dim_velocity, red_dim_pressure, B_velocity, B_pressure)

rhs_reduced = zeros(red_dim_velocity+red_dim_pressure,1);
rhs_reduced(1:red_dim_velocity) = B_velocity' * rhs(1:params.ndofs);
rhs_reduced(red_dim_velocity+1:red_dim_velocity+red_dim_pressure) = ...
    B_pressure' * rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
stifness_matrix_reduced = zeros(red_dim_velocity+red_dim_pressure);
stifness_matrix_reduced(1:red_dim_velocity,1:red_dim_velocity) = ...
    B_velocity' * stifness_matrix(1:params.ndofs,1:params.ndofs) * B_velocity;
stifness_matrix_reduced(1:red_dim_velocity, ...
    red_dim_velocity + 1 : red_dim_velocity + red_dim_pressure) ...
    = B_velocity' * stifness_matrix(1:params.ndofs,params.ndofs + 1 : ...
    params.ndofs + paramsP.ndofs) * B_pressure;
stifness_matrix_reduced(red_dim_velocity+1:red_dim_velocity + ...
    red_dim_pressure, 1:red_dim_velocity) = B_pressure' * ...
    stifness_matrix(params.ndofs + 1:params.ndofs + paramsP.ndofs, ...
    1:params.ndofs) * B_velocity;
stifness_matrix_reduced = sparse(stifness_matrix_reduced);
params_reduced.ndofs = red_dim_velocity;
paramsP_reduced.ndofs = red_dim_pressure;
params_reduced.dofs = zeros(red_dim_velocity,1);
paramsP_reduced.dofs = zeros(red_dim_pressure,1);
params_reduced.bilinear_side = ...
    stifness_matrix_reduced(1:red_dim_velocity,1:red_dim_velocity);
params_reduced.bilinear_side_pressure_terms = stifness_matrix_reduced...
    (1:red_dim_velocity, red_dim_velocity+1:red_dim_velocity+red_dim_pressure);
params_reduced.linear_side = rhs_reduced(1:red_dim_velocity);
params_reduced.rhs_continuity = rhs_reduced(red_dim_velocity+1:...
    red_dim_velocity+red_dim_pressure);

params_reduced.dimrange = params.dimrange;
paramsP_reduced.dimrange = paramsP.dimrange;
params_reduced.pdeg = params.pdeg;
paramsP_reduced.pdeg = paramsP.pdeg;
params_reduced.grid = grid;
paramsP_reduced.grid = grid;
params_reduced.ndofs_per_element = red_dim_velocity;
paramsP_reduced.ndofs_per_element = red_dim_pressure;

[ params_reduced, paramsP_reduced] = solve_plot_solution_schur...
    ( params_reduced, paramsP_reduced, grid, rhs_reduced, ...
    stifness_matrix_reduced);

% dofs = minres(stifness_matrix_reduced,rhs_reduced,1e-11,1e4);
% 
% params_reduced.dofs = dofs(1:red_dim_velocity);
% paramsP_reduced.dofs = dofs(red_dim_velocity+1:red_dim_velocity+red_dim_pressure);

end