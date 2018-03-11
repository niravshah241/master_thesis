% grid making and saving
% pdegrid
% save('mygrid','p','e','t')
clc
%% Parameter generation
disp('Generating reference parameter set')
para1 = 1; %viscocity
para2 = 0; %dirichlet value
parameter_reference_set = [para1 para2];
params.reference_parameter = parameter_reference_set;
disp('Parameter reference set created')

%% Parameter generation
disp('Generating training parameter set')
para1 = [1 1.1 3]; %viscocity
para2 = [0 0.1 3]; %dirichlet value
parameter_training_set = gen_parameters( para1,para2);
disp('Parameter training set generation finished')

%% Actual grid
% ACTUAL GRID DO NOT DELETE

params.mesh_number = 3;
params.gridtype = 'triagrid';
params.grid_initfile = ['mygridnirav', num2str(params.mesh_number), '.mat'];
% params.bnd_rect_corner1=[-1,-1;-eps,eps]'; % for analytical
% params.bnd_rect_corner2=[eps,1+eps;eps,1-3*10^14*eps]';% for analytical ex.
% params.bnd_rect_corner1=[-1,-1;100,10]'; % for benchmark problem
% params.bnd_rect_corner2=[2,2;100,10-eps]'; % for benchmark problem
params.bnd_rect_corner1=[-1,-1;1-eps,3*10^14*eps]'; % for standard
params.bnd_rect_corner2=[eps,1+eps;1+eps,1-eps]'; % for standard
params.bnd_rect_index=[-1,-2];
grid=construct_grid(params);
show_sparsity = false; % Bool variable which plots sparsity pattern of
% % assembled matrix is set to true else(i.e. false) the sparsity pattern is not shown
params.show_sparsity = show_sparsity;
paramsP.show_sparsity = show_sparsity;

%ACTUAL GRID OVER

%% Test grid

%ONLY FOR TEST GRID

% params.xrange = [0,1];
% params.yrange = [0,1];
% params.xnumintervals = 10;
% params.ynumintervals = 10;
% params.bnd_rect_corner1=[-1,-1;-eps,eps]'; % for analytical
% params.bnd_rect_corner2=[2,2;eps,1-0.06]';% for analytical ex.
% % params.bnd_rect_corner1=[-1,-1;1-eps,0+3*10^14*eps]';
% % params.bnd_rect_corner2=[eps,1+eps;1+eps,1-eps]';
% params.bnd_rect_index=[-1,-2];
% params.gridtype = 'triagrid';
% grid = construct_grid(params);
% show_sparsity = false; % Bool variable which plots sparsity pattern of
% %assembled matrix is set to true else(i.e. false) the sparsity pattern is not shown
% params.show_sparsity = show_sparsity;
% paramsP.show_sparsity = show_sparsity;

%TEST GRID OVER

%% Plotting of grid
disp('Please check the grid')
figure()
plot(grid);
title('Grid')
% pause();
% close all

%% Variables setting
params.pdeg = 2;
paramsP.pdeg = params.pdeg-1;%taylor hood element
params.dimrange = 2;
paramsP.dimrange = 1;
params.grid = grid;
paramsP.grid = grid;

nrep = [3 6 10 15];

params.ndofs_per_element = nrep(params.pdeg) * params.dimrange;
params.ndofs = params.ndofs_per_element * grid.nelements;
params.dofs = zeros(params.ndofs,1);

paramsP.ndofs_per_element = nrep(paramsP.pdeg) * paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element * grid.nelements;
paramsP.dofs = zeros(paramsP.ndofs,1);

df_info = ldginfo(params,grid);
df = ldgdiscfunc(df_info);
display(df);
df_infoP = ldginfo(paramsP,grid);
df = ldgdiscfunc(df_info);
display(df);
qdeg = 3;

params.kinematic_viscosity = @(params) params.reference_parameter(1)*1e-6;
mu = params.kinematic_viscosity(params);
c11 = 1e1*mu*1e6;% penalty parameter, must be large enough for coercivity

%% Assembly of stiffness matrix

disp('Solving reference problem')
disp('Assembling stifness matrix')

tic
[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 );
time_matrix_assembly = toc;
disp(['Time taken for assembling stifness matrix ',num2str(time_matrix_assembly)])
stifness_matrix_reference = stifness_matrix;
rhs_reference = rhs;

%% Snapshot generation
disp('Entering in snapshot generation')
params.snapshots_matrix = zeros(params.ndofs,size(parameter_training_set,1));
paramsP.snapshots_matrix = zeros(paramsP.ndofs,size(parameter_training_set,1));
params_reference = params;
paramsP_reference = paramsP;

for i = 1:1:size(parameter_training_set,1)
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_training_set(i,:);
    params.kinematic_viscosity = @(params) params.parameter_training_set(1)*1e-6;
    mu = params.kinematic_viscosity(params);
    c11 = 1e4; %1e1*mu*1e6;% penalty parameter, must be large enough for coercivity
    
    %% AFFINE assembly of stiffness matrix
    
    %[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    %    ( params, paramsP, grid, qdeg, mu, c11 );
    
    theta_1 = params.parameter_training_set(1) - params.reference_parameter(1);
    %viscocity
    theta_2 = params.parameter_training_set(2) / params.reference_parameter(2);
    %dirichlet value
    
    disp('Assembling stifness matrix affine')
    tic
    stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
        stifness_matrix(1:params.ndofs,1:params.ndofs) + theta_1 * ...
        (params_reference.bilinear_res1 - params_reference.bilinear_res2 ...
        - params_reference.bilinear_res2');
    
    rhs(1:params.ndofs) = rhs(1:params.ndofs) - theta_1 * ...
        params_reference.linear_res4 + (theta_2 - 1) * ...
        (params_reference.linear_res3 - params_reference.linear_res4);
    
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        theta_2*rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    
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

%% Proper Orthogonal Decomposition
n_s = size(params.snapshots_matrix,2); % number of snapshots
% params.snapshots_matrix = rand(params.ndofs,n_s);
red_dim_velocity = 15;
red_dim_pressure = 6;
min_eigen = 1e-10;
params.qdeg = qdeg;
paramsP.qdeg = qdeg;

disp('Entering in pod')

[ pod_res_params, pod_res_paramsP, B, B_velocity, B_pressure, ...
    stifness_matrix_reduced, params_reduced, paramsP_reduced] = ...
    pod( params, paramsP, grid, red_dim_velocity, red_dim_pressure, ...
    min_eigen, stifness_matrix, rhs);

%% Testing
disp('Generating test parameter set')
para_test_1 = [1 10 3];
para_test_2 = [0 0.1 2];
parameter_test_set = gen_test_parameters(para_test_1,para_test_2);
disp('Test parameter set generated')
error_velocity = zeros(size(parameter_test_set,1),1);
error_pressure = zeros(size(parameter_test_set,1),1);

error_velocity_vector = zeros(size(parameter_test_set,1),1);
error_pressure_vector = zeros(size(parameter_test_set,1),1);

disp('Entering in error calculation')
for i = 1:1:size(parameter_test_set,1)
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_test_set(i,:);
    params.kinematic_viscosity = @(params) params.parameter_training_set(1)*1e-6;
    mu = params.kinematic_viscosity(params);
    c11 = 1e1*mu*1e6;% penalty parameter, must be large enough for coercivity
    
    %% Assembly of stiffness matrix
    
%     [ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
%         ( params, paramsP, grid, qdeg, mu, c11 );
    
    %% Stokes problem
    
    
    theta_1 = params.parameter_training_set(1) - params.reference_parameter(1);
    theta_2 = params.parameter_training_set(2) / params.reference_parameter(2);
    
    disp('Assembling stifness matrix affine')
    tic
    stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
        stifness_matrix(1:params.ndofs,1:params.ndofs) + theta_1 * ...
        (params_reference.bilinear_res1 - params_reference.bilinear_res2 - ...
        params_reference.bilinear_res2');
    
    rhs(1:params.ndofs) = rhs(1:params.ndofs) - theta_1 * ...
        params_reference.linear_res4 + (theta_2 - 1) * (params_reference.linear_res3 ...
        - params_reference.linear_res4);
    
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        theta_2*rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    
    time_matrix_assembly_affine = toc;
    disp(['Time taken for affine assembling stifness matrix ',...
        num2str(time_matrix_assembly_affine)])
    
    tic;
    [ params, paramsP, achieved_residual_tol_schur] =...
        solve_plot_solution_schur( params, paramsP, grid, rhs, stifness_matrix);
    time_schur = toc;
    
    error_velocity_vector(i) = error_velocity_rbasis( params, params_reduced, B_velocity, ...
        grid, qdeg);
    error_pressure_vector(i) = error_velocity_rbasis( paramsP, paramsP_reduced, B_pressure, ...
        grid, qdeg);
    close all
end

% Error in energy norm

error_energy_velocity_reduced = B_velocity'*stifness_matrix(1:params.ndofs,...
    1:params.ndofs)*B_velocity*params_reduced.dofs + ...
    B_velocity'*stifness_matrix(1:params.ndofs,params.ndofs + ...
    1:params.ndofs+paramsP.ndofs)*B_pressure*paramsP_reduced.dofs;
error_energy_pressure_reduced = B_pressure'*...
    stifness_matrix(params.ndofs+1:params.ndofs + ...
    paramsP.ndofs,1:params.ndofs)*B_velocity*params_reduced.dofs;
error_energy_velocity = B_velocity' * stifness_matrix(1:params.ndofs,1:params.ndofs) * ...
    params.dofs + B_velocity' * stifness_matrix(1:params.ndofs,params.ndofs + ...
    1:params.ndofs+paramsP.ndofs) * paramsP.dofs - error_energy_velocity_reduced;
error_energy_velocity = norm(error_energy_velocity,2);
error_energy_pressure = B_pressure' * stifness_matrix(params.ndofs+1:params.ndofs ...
    + paramsP.ndofs,1:params.ndofs)*params.dofs - error_energy_pressure_reduced;
error_energy_pressure = norm(error_energy_pressure,2);