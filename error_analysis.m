function [error_l2, error_energy, velocity_error, pressure_error] = error_analysis(params, paramsP, ...
    stifness_matrix_reference, rhs_reference, parameter_test_set, B_velocity, ...
    B_pressure, grid, para_test_1, para_test_2, linear_side_reference, ...
    red_dim_velocity, red_dim_pressure, reference_factor)

error_energy = zeros(size(parameter_test_set,1),1);
error_l2 = zeros(size(parameter_test_set,1),1);

for i = 1:1:size(parameter_test_set,1)
    disp(['Test parameter number ',num2str(i),' of ', ...
        num2str(para_test_1(3)*para_test_2(3))])
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_test_set(i,:);
    mu = params.parameter_training_set(1);
    
    %% Stokes problem
    
    
    theta_1 = params.parameter_training_set(1) / params.reference_parameter(1);
    %viscocity
    theta_2 = params.parameter_training_set(2) / params.reference_parameter(2);
    %dirichlet value
    
    disp('Assembling stifness matrix affine')
    tic
    stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
        theta_1 * stifness_matrix(1:params.ndofs,1:params.ndofs);
    
    rhs(1:params.ndofs) = rhs(1:params.ndofs) + ...
        (params.parameter_training_set(1) * theta_2 - ...
        params.reference_parameter(1)) * (reference_factor * ...
        linear_side_reference.term_3 - linear_side_reference.term_4);
    
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        theta_2 * rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    
    time_matrix_assembly_affine = toc;
    disp(['Time taken for affine assembling stifness matrix ',...
        num2str(time_matrix_assembly_affine)])
    
    tic;
    [ params, paramsP, achieved_residual_tol_schur] =...
        solve_plot_solution_schur( params, paramsP, grid, rhs, stifness_matrix);
    time_schur = toc;
    
    [params_reduced, paramsP_reduced, stifness_matrix_reduced, rhs_reduced] = ...
        galerkin_formulation(stifness_matrix, rhs, params, paramsP, grid, ...
        red_dim_velocity, red_dim_pressure, B_velocity, B_pressure);
    
    params_reduced.dofs = B_velocity * params_reduced.dofs;
    velocity_error = params.dofs - params_reduced.dofs;
    paramsP_reduced.dofs = B_pressure * paramsP_reduced.dofs;
    pressure_error = paramsP.dofs - paramsP_reduced.dofs;
    
    error_l2(i) = velocity_error' * ldg_mass_matrix(params,grid,params)...
        *velocity_error + pressure_error' * ldg_mass_matrix...
        (paramsP,grid,paramsP)*pressure_error;
    dofs = [params.dofs;paramsP.dofs];
    error_l2(i) = error_l2(i) / (params.dofs' * ldg_mass_matrix...
        (params,grid,params) * params.dofs + paramsP.dofs' * ...
        ldg_mass_matrix(paramsP,grid,paramsP) * paramsP.dofs);
    error_energy(i) = [velocity_error;pressure_error]' * ...
        stifness_matrix * [velocity_error;pressure_error]; 
    error_energy(i) = error_energy(i) / (dofs'*stifness_matrix * dofs);
    close all
end

end