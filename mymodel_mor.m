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
para1 = [1 3 5]; %viscocity
para2 = [1 2 5]; %dirichlet value
parameter_training_set = gen_parameters( para1,para2);
disp('Parameter training set generation finished')

%% Snapshot generation
disp('Entering in snapshot generation')
params.snapshots_matrix = zeros(params.ndofs,size(parameter_training_set,1));
paramsP.snapshots_matrix = zeros(paramsP.ndofs,size(parameter_training_set,1));
params_reference = params;
paramsP_reference = paramsP;

[params, paramsP] = snapshot_generation( params, paramsP, para1, ...
    para2, stifness_matrix_reference, rhs_reference, grid, ...
    parameter_training_set, linear_side_reference);

%% Proper Orthogonal Decomposition
n_s = size(params.snapshots_matrix,2); % number of snapshots
% params.snapshots_matrix = rand(params.ndofs,n_s);
red_dim_velocity = 20;
red_dim_pressure = 20;
params.qdeg = qdeg;
paramsP.qdeg = qdeg;
min_eigen_velocity = 1e-12;
min_eigen_pressure = 1e-12;

disp('Entering in pod')

[ pod_res_params, pod_res_paramsP, B_velocity, B_pressure, ...
    red_dim_velocity, red_dim_pressure] = pod( params, paramsP, grid, ...
    red_dim_velocity, red_dim_pressure, min_eigen_pressure, ...
    min_eigen_velocity, stifness_matrix, rhs);

%% Testing
disp('Generating test parameter set')
para_test_1 = [1 2 10];
para_test_2 = [1 2 10];
parameter_test_set = gen_test_parameters(para_test_1,para_test_2);
disp('Test parameter set generated')

disp('Entering in error calculation')

[error_l2, error_energy] = error_analysis(params, paramsP, ...
    stifness_matrix_reference, rhs_reference, parameter_test_set, B_velocity, ...
    B_pressure, grid, para_test_1, para_test_2, linear_side_reference, ...
    red_dim_velocity, red_dim_pressure);

figure()
[xq,yq] = meshgrid(para_test_1(1):0.1:para_test_1(2),...
    para_test_2(1):0.1:para_test_2(2));
vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
    error_l2,xq,yq);
mesh(xq,yq,vq)
hold on
plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_l2,'o')
axis tight
xlabel('Viscocity')
ylabel('Dirichlet value')
zlabel('Error L^2')
title('Error L^2 over parameter space')

figure()
[xq,yq] = meshgrid(para_test_1(1):0.1:para_test_1(2),...
    para_test_2(1):0.1:para_test_2(2));
vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
    error_energy,xq,yq);
mesh(xq,yq,vq)
hold on
plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_energy,'o')
axis tight
xlabel('Viscocity')
ylabel('Dirichlet value')
zlabel('Error energy')
title('Error energy over parameter space')

% online phase

parameter_online = [2 1.5];
params.parameter_online = parameter_online;

tic
[params, paramsP, params_reduced, paramsP_reduced] = online_phase...
    (params, paramsP, grid, stifness_matrix_reference, rhs_reference, ...
    B_velocity, B_pressure, linear_side_reference, red_dim_velocity, ...
    red_dim_pressure);
online_time = toc;
disp(['Time for online phase: ',num2str(online_time)])

tic
[params, paramsP] = dg_solution( params, paramsP, grid, rhs_reference, ...
    stifness_matrix_reference, linear_side_reference);
offline_time = toc;
disp(['Time for ofline phase: ',num2str(offline_time)])