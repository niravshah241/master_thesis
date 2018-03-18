% grid making and saving
% pdegrid
% save('mygrid','p','e','t')
clc
%% Parameter-reference generation
disp('Generating reference parameter set')
para1 = 1; %viscocity
para2 = 1; %dirichlet value
parameter_reference_set = [para1 para2];
params.reference_parameter = parameter_reference_set;
disp('Parameter reference set created')

%% Actual grid
% ACTUAL GRID DO NOT DELETE

params.mesh_number = 1;
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
dfP = ldgdiscfunc(df_infoP);
display(dfP);
qdeg = 3;

params.kinematic_viscosity = @(params) params.reference_parameter(1)*1e-6;
mu = params.kinematic_viscosity(params);
reference_factor = 1e6;
c11 = reference_factor*mu;% penalty parameter, must be large enough for coercivity

%% Assembly of stiffness matrix

disp('Solving reference problem')
disp('Assembling stifness matrix')

tic
[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 );
linear_side_reference.term_1 = params.linear_res1;
linear_side_reference.term_2 = params.linear_res2;
linear_side_reference.term_3 = params.linear_res3;
linear_side_reference.term_4 = params.linear_res4;
time_matrix_assembly = toc;
disp(['Time taken for assembling stifness matrix ',num2str(time_matrix_assembly)])
stifness_matrix_reference = stifness_matrix;
rhs_reference = rhs;

%% Parameter-training generation
disp('Generating training parameter set')
para1 = [1e-10 1e5 3]; %viscocity
para2 = [1 1e8 3]; %dirichlet value
parameter_training_set = gen_parameters( para1,para2);
disp('Parameter training set generation finished')

%% Snapshot generation
disp('Entering in snapshot generation')
params.snapshots_matrix = zeros(params.ndofs,size(parameter_training_set,1));
paramsP.snapshots_matrix = zeros(paramsP.ndofs,size(parameter_training_set,1));
params_reference = params;
paramsP_reference = paramsP;

for i = 1:1:size(parameter_training_set,1)
    disp(['Training parameter number ',num2str(i),' of ', ...
        num2str(para1(3)*para2(3))])
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_training_set(i,:);
    params.kinematic_viscosity = @(params) params.parameter_training_set(1)*1e-6;
    mu = params.kinematic_viscosity(params);
    %c11 = 1e1*mu*1e3;% penalty parameter, must be large enough for coercivity
    
    %% AFFINE assembly of stiffness matrix
    
    %[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    %    ( params, paramsP, grid, qdeg, mu, c11 );
    
    theta_1 = params.parameter_training_set(1) / params.reference_parameter(1);
    %viscocity
    theta_2 = params.parameter_training_set(2) / params.reference_parameter(2);
    %dirichlet value
    
    disp('Assembling stifness matrix affine')
    tic
    stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
        theta_1 * stifness_matrix(1:params.ndofs,1:params.ndofs);
    
    rhs(1:params.ndofs) = rhs(1:params.ndofs) + (theta_1 * theta_2 - 1)...
        * (linear_side_reference.term_3 - linear_side_reference.term_4);
    
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
red_dim_velocity = 2 * grid.nelements;
red_dim_pressure = 3 * grid.nelements;
params.qdeg = qdeg;
paramsP.qdeg = qdeg;
min_eigen_velocity = 0;
min_eigen_pressure = 0;

disp('Entering in pod')


[ pod_res_params, pod_res_paramsP, B_velocity, B_pressure, ...
    red_dim_velocity, red_dim_pressure] = pod( params, paramsP, grid, ...
    red_dim_velocity, red_dim_pressure, min_eigen_pressure, ...
    min_eigen_velocity, stifness_matrix, rhs);

%% Testing
disp('Generating test parameter set')
para_test_1 = [1e-10 1e5 3];
para_test_2 = [1 1e8 3];
parameter_test_set = gen_test_parameters(para_test_1,para_test_2);
disp('Test parameter set generated')
error_velocity = zeros(size(parameter_test_set,1),1);
error_pressure = zeros(size(parameter_test_set,1),1);

error_velocity_vector = zeros(size(parameter_test_set,1),1);
error_pressure_vector = zeros(size(parameter_test_set,1),1);
error_energy_velocity = ones(size(parameter_test_set,1),1);
error_energy_pressure = ones(size(parameter_test_set,1),1);
error_energy = zeros(size(parameter_test_set,1),1);

disp('Entering in error calculation')
for i = 1:1:size(parameter_test_set,1)
    disp(['Test parameter number ',num2str(i),' of ', ...
        num2str(para_test_1(3)*para_test_2(3))])
    stifness_matrix = stifness_matrix_reference;
    rhs = rhs_reference;
    params.parameter_training_set = parameter_test_set(i,:);
    params.kinematic_viscosity = @(params) params.parameter_training_set(1)*1e-6;
    mu = params.kinematic_viscosity(params);
    %c11 = 1e3*mu*1e3;% penalty parameter, must be large enough for coercivity
    
    %% Assembly of stiffness matrix
    
    %     [ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    %         ( params, paramsP, grid, qdeg, mu, c11 );
    
    %% Stokes problem
    
    
    theta_1 = params.parameter_training_set(1) / params.reference_parameter(1);
    %viscocity
    theta_2 = params.parameter_training_set(2) / params.reference_parameter(2);
    %dirichlet value
    
    disp('Assembling stifness matrix affine')
    tic
    stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
        theta_1 * stifness_matrix(1:params.ndofs,1:params.ndofs);
    
    rhs(1:params.ndofs) = rhs(1:params.ndofs) + (theta_1 * theta_2 - 1)...
        * (linear_side_reference.term_3 - linear_side_reference.term_4);
    
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        theta_2*rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    
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
    
    error_velocity_vector(i) = error_velocity_rbasis( params, params_reduced, B_velocity, ...
        grid, qdeg);
    error_pressure_vector(i) = error_velocity_rbasis( paramsP, paramsP_reduced, B_pressure, ...
        grid, qdeg);
    
    % Error in energy norm
    
    error_energy(i) = [params.dofs',paramsP.dofs'] * stifness_matrix * ...
        [params.dofs;paramsP.dofs] - [params_reduced.dofs',paramsP_reduced.dofs'] * ...
        stifness_matrix_reduced * [params_reduced.dofs;paramsP_reduced.dofs];
    
    close all
    
end

params_reduced.dofs
paramsP_reduced.dofs

% tria_index = 1;
% glob = [0 1];
% tria_index = 1;
% [res_velocity] = online_solve(glob, tria_index, params, params_reduced, ...
%     grid, B_velocity);
% [res_pressure] = online_solve(glob, tria_index, paramsP, paramsP_reduced, ...
%     grid, B_pressure);