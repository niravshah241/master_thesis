function [params, paramsP] = snapshot_generation( params, paramsP, para1, ...
    para2, stifness_matrix_reference, rhs_reference, grid, ...
    parameter_training_set, linear_side_reference, reference_factor);

for i = 1:1:size(parameter_training_set,1)
    disp(['Training parameter number ',num2str(i),' of ', ...
        num2str(para1(3)*para2(3))])
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_training_set(i,:);
    mu = params.parameter_training_set(1);
    %c11 = 1e1*mu*1e3;% penalty parameter, must be large enough for coercivity
    
    %% AFFINE assembly of stiffness matrix
    
    %[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    %    ( params, paramsP, grid, qdeg, mu, c11 );
    
    theta_1 = params.parameter_training_set(1) / ...
        params.reference_parameter(1);
    %viscocity
    theta_2 = params.parameter_training_set(2) / ...
        params.reference_parameter(2);
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
    
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        theta_2 * rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    
    time_matrix_assembly_affine = toc;
    disp(['Time taken for affine assembling stifness matrix ',...
        num2str(time_matrix_assembly_affine)])
    
    %% Stokes problem
    
    tic;
    [ params, paramsP, achieved_residual_tol_schur] =...
        solve_plot_solution_schur( params, paramsP, grid, rhs, stifness_matrix);
    time_schur = toc;
    
    required_residual_tol = 0;%achieved_residual_tol_schur; % allowable residual
    max_iter = 1e5; % maximum number of iterations
    
    % tic;
    % [ params, paramsP, flag, achieved_residual_tol, actual_iter] = solve_plot_solution...
    %     ( params, paramsP, grid, rhs, stifness_matrix, required_residual_tol, max_iter);
    % times_solver = toc;
    
    params.snapshots_matrix(:,i) = params.dofs;
    paramsP.snapshots_matrix(:,i) = paramsP.dofs;
    close all
end
end