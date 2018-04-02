clc
clear all
%% Parameter-reference generation
disp('Generating reference parameter set')
para1 = 1; %viscocity
para2 = 1; %dirichlet value
parameter_reference_set = [para1 para2];
params.reference_parameter = parameter_reference_set;
disp('Parameter reference set created')


%% Parameter-training generation
disp('Generating training parameter set')
para1 = [1 2 10]; %viscocity
para2 = [10e-2 20e-2 10]; %dirichlet value
parameter_training_set = gen_parameters( para1,para2);
disp('Parameter training set generation finished')

%% Testing
disp('Generating test parameter set')
para_test_1 = [1 2 5];
para_test_2 = [10e-2 20e-2 5];
parameter_test_set = gen_test_parameters(para_test_1,para_test_2);
disp('Test parameter set generated')

%% Online parameter

parameter_online = [1 15e-2];
params.parameter_online = parameter_online;

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

mu = params.reference_parameter(1);
reference_factor = 1e6;
c11 = reference_factor*mu;% penalty parameter, must be large enough for coercivity

%% Assembly of stiffness matrix

disp('Solving reference problem')
disp('Assembling stifness matrix')

tic
[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 );
time_matrix_assembly = toc;
disp(['Time taken for assembling stifness matrix ',num2str(time_matrix_assembly)])
linear_side_reference.term_1 = params.linear_res1;
linear_side_reference.term_2 = params.linear_res2;
linear_side_reference.term_3 = params.linear_res3;
linear_side_reference.term_4 = params.linear_res4;
stifness_matrix_reference = stifness_matrix;
rhs_reference = rhs;


%% Snapshot generation
disp('Entering in snapshot generation')
params.snapshots_matrix = zeros(params.ndofs,size(parameter_training_set,1));
paramsP.snapshots_matrix = zeros(paramsP.ndofs,size(parameter_training_set,1));
params_reference = params;
paramsP_reference = paramsP;

[params, paramsP] = snapshot_generation( params, paramsP, para1, ...
    para2, stifness_matrix_reference, rhs_reference, grid, ...
    parameter_training_set, linear_side_reference, reference_factor);

x = 1:5:para1(3)*para2(3);

error_l2_average = zeros(length(x),1);
error_energy_average = zeros(length(x),1);
online_time_vector = zeros(length(x),1);

for i = 1:1:length(x)
    %% Proper Orthogonal Decomposition
    n_s = size(params.snapshots_matrix,2); % number of snapshots
    % params.snapshots_matrix = rand(params.ndofs,n_s);
    red_dim_velocity = x(i);
    red_dim_pressure = x(i);
    params.qdeg = qdeg;
    paramsP.qdeg = qdeg;
    min_eigen_velocity = 1e-12;
    min_eigen_pressure = 1e-12;
    
    disp(['Entering in pod ',num2str(i),' of ',num2str(length(x))])
    
    [ pod_res_params, pod_res_paramsP, B_velocity, B_pressure, ...
        red_dim_velocity, red_dim_pressure] = pod( params, paramsP, grid, ...
        red_dim_velocity, red_dim_pressure, min_eigen_pressure, ...
        min_eigen_velocity, stifness_matrix, rhs);
    
    %% Error analysis
    
    disp('Entering in error calculation')
    
    [error_l2, error_energy, velocity_error, pressure_error] = ...
        error_analysis(params, paramsP, stifness_matrix_reference, ...
        rhs_reference, parameter_test_set, B_velocity, B_pressure, ...
        grid, para_test_1, para_test_2, linear_side_reference, ...
        red_dim_velocity, red_dim_pressure, reference_factor);
    
    figure()
    [xq,yq] = meshgrid(parameter_training_set(:,1),...
        parameter_training_set(:,2));
    vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
        error_l2,xq,yq);
    mesh(xq,yq,vq)
    hold on
    plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_l2,'o')
    axis tight
    xlabel('Viscocity')
    ylabel('Dirichlet value')
    zlabel('Relative error L^2')
    title('Error L^2 over parameter space')
    
    figure()
    [xq,yq] = meshgrid(parameter_training_set(:,1),...
        parameter_training_set(:,2));
    vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
        error_energy,xq,yq);
    mesh(xq,yq,vq)
    hold on
    plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_energy,'o')
    axis tight
    xlabel('Viscocity')
    ylabel('Dirichlet value')
    zlabel('Relative error energy')
    title('Error energy over parameter space')
    error_l2_average(i) = mean(error_l2);
    error_energy_average(i) = mean(error_energy);
    tic
    [params, paramsP, params_reduced, paramsP_reduced] = online_phase...
        (params, paramsP, grid, stifness_matrix_reference, rhs_reference, ...
        B_velocity, B_pressure, linear_side_reference, red_dim_velocity, ...
        red_dim_pressure, reference_factor);
    online_time_vector(i) = toc;
    disp(['Time for online phase: ',num2str(online_time_vector(i))])
end

figure()
plot(x,error_l2_average);
title('L^2 error')
xlabel('Size of reduced basis')
ylabel('Relative L^2 error')
axis tight

figure()
plot(x,error_energy_average);
title('Relative energy error')
xlabel('Size of reduced basis')
ylabel('energy error')
axis tight

%% Test for affine stifness matrix vs actual stifness matrix

mu = params.parameter_online(1);
c11 = reference_factor * mu;
params.parameter_training_set = params.parameter_online;

[ params, paramsP, rhs_assembled, stifness_matrix_assembled] = ...
    assemble_stifness_matrix( params, paramsP, grid, qdeg, mu, c11 );

stifness_matrix = stifness_matrix_reference;
rhs = rhs_reference;

theta_1 = params.parameter_online(1) / params.reference_parameter(1);
%viscocity
theta_2 = params.parameter_online(2) / params.reference_parameter(2);
%dirichlet value

disp('Assembling stifness matrix affine')

stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
    theta_1 * stifness_matrix(1:params.ndofs,1:params.ndofs);

mu_online = params.parameter_online(1);
mu_reference = params.reference_parameter(1);

rhs(1:params.ndofs) = rhs(1:params.ndofs) + ...
    (mu_online * theta_2 - mu_reference) * (reference_factor * ...
    linear_side_reference.term_3 - linear_side_reference.term_4);

rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
    theta_2 * rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);

stifness_matrix_affine = stifness_matrix;
rhs_affine = rhs;
error_stifness_matrix_max = max(max(stifness_matrix_assembled - ...
    stifness_matrix_affine));
error_stifness_matrix_2norm = norm(full(stifness_matrix_assembled) - ...
    full(stifness_matrix_affine),2);
error_rhs_max = max(abs(rhs_assembled - rhs_affine));
error_rhs_2norm = norm(full(rhs_assembled) - full(rhs_affine),2);

%% Online phase

tic
[params, paramsP, params_reduced, paramsP_reduced] = online_phase...
    (params, paramsP, grid, stifness_matrix_reference, rhs_reference, ...
    B_velocity, B_pressure, linear_side_reference, red_dim_velocity, ...
    red_dim_pressure, reference_factor);
online_time = toc;
disp(['Time for online phase: ',num2str(online_time)])

error_velocity.dofs = params.dofs - params_reduced.dofs;
error_velocity.pdeg = params.pdeg;
error_velocity.dimrange = params.dimrange;
error_velocity.ndofs_per_element = params.ndofs_per_element;
error_velocity.ndofs = params.ndofs;
error_velocity.grid = grid;
error_pressure.dimrange = paramsP.dimrange;
error_pressure.ndofs_per_element = paramsP.ndofs_per_element;
error_pressure.ndofs = paramsP.ndofs;
error_pressure.dofs = paramsP.dofs - paramsP_reduced.dofs;
error_pressure.pdeg = paramsP.pdeg;
error_pressure.grid = grid;
for i=1:1:error_velocity.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(error_velocity,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting',num2str(i),' degree of freedom (error velocity)'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Error velocity in x direction'])
    else
        title(['Error velocity in y direction'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,grid,error_velocity);
    plot(grid);
end

for i=1:1:error_pressure.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(error_pressure,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (Error pressure)'])
    %subplot(paramsP.dimrange,1,i)
    title('Error in pressure')
    axis equal
    axis tight
    ldg_plot(sdf,grid,error_pressure);
    plot(grid);
end

tic
[params, paramsP] = dg_solution( params, paramsP, grid, ...
    reference_factor, qdeg);
offline_time = toc;
disp(['Time for offline phase: ',num2str(offline_time)])

figure()
plot(x,online_time_vector/offline_time,'*')
axis tight
title('Time saving')
xlabel('Size of reduced basis space')
ylabel('Time for online computation/Time for offline computation')

figure()
[xq,yq] = meshgrid(parameter_training_set(:,1),...
    parameter_training_set(:,2));
vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
    error_l2,xq,yq);
mesh(xq,yq,vq)
hold on
plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_l2,'o')
axis tight
xlabel('Viscocity')
ylabel('Dirichlet value')
zlabel('Relative error L^2')
title('Error L^2 over parameter space')

figure()
[xq,yq] = meshgrid(parameter_training_set(:,1),...
    parameter_training_set(:,2));
vq = griddata(parameter_test_set(:,1),parameter_test_set(:,2),...
    error_energy,xq,yq);
mesh(xq,yq,vq)
hold on
plot3(parameter_test_set(:,1),parameter_test_set(:,2),error_energy,'o')
axis tight
xlabel('Viscocity')
ylabel('Dirichlet value')
zlabel('Relative error energy')
title('Error energy over parameter space')
error_l2_average(i) = mean(error_l2);
error_energy_average(i) = mean(error_energy);
tic
[params, paramsP, params_reduced, paramsP_reduced] = online_phase...
    (params, paramsP, grid, stifness_matrix_reference, rhs_reference, ...
    B_velocity, B_pressure, linear_side_reference, red_dim_velocity, ...
    red_dim_pressure, reference_factor);
online_time_vector(i) = toc;
disp(['Time for online phase: ',num2str(online_time_vector(i))])