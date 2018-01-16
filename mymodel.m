%grid making and saving
%pdegrid
%save('mygrid','p','e','t')
clc
%ACTUAL GRID DO NOT DELETE

params.mesh_number = 1;
model.gridtype = 'triagrid';
model.grid_initfile = ['mygridnirav', num2str(params.mesh_number), '.mat'];
model.bnd_rect_corner1=[-1,-1;1-eps,0+3*10^14*eps]';
model.bnd_rect_corner2=[2,2;1+eps,1-eps]';
model.bnd_rect_corner2=[eps,1+eps;1+eps,1-eps]';
model.bnd_rect_index=[-1,-2];
grid=construct_grid(model);
show_sparsity = false; % Bool variable which plots sparsity pattern of 
% assembled matrix is set to true else(i.e. false) the sparsity pattern is not shown
params.show_sparsity = show_sparsity;
paramsP.show_sparsity = show_sparsity;

%ACTUAL GRID OVER

%ONLY FOR TEST GRID

% params.xrange = [0,1];
% params.yrange = [0,1];
% params.xnumintervals = 3;
% params.ynumintervals = 3;
% params.bnd_rect_corner1=[-1,-1;1-eps,1*10^15*eps]';
% params.bnd_rect_corner2=[2,2;1+eps,1-eps]';
% %params.bnd_rect_corner2=[1-eps,eps;1+eps,1-eps]';
% params.bnd_rect_index=[-1,-2];
% grid = triagrid(params);
% show_sparsity = false; % Bool variable which plots sparsity pattern of 
% %assembled matrix is set to true else(i.e. false) the sparsity pattern is not shown
% params.show_sparsity = show_sparsity;
% paramsP.show_sparsity = show_sparsity;

%TEST GRID OVER

disp('Please check the grid')
figure()
plot(grid);
title('Grid')
%pause();
%close all
params.pdeg = 1;
paramsP.pdeg = 1;%params.pdeg-1;%taylor hood element
params.dimrange = 2;
paramsP.dimrange = 1;
params.grid = grid;
paramsP.grid = grid;

nrep=[3 6 10 15];

params.ndofs_per_element= nrep(params.pdeg)*params.dimrange;
params.ndofs = params.ndofs_per_element*grid.nelements;
params.dofs = ones(params.ndofs,1);

paramsP.ndofs_per_element= nrep(paramsP.pdeg)*paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element*grid.nelements;
paramsP.dofs = ones(paramsP.ndofs,1);

df_info=ldginfo(params,grid);
df=ldgdiscfunc(df_info);
display(df);
qdeg=3;
params.mu=4;
params.kinematic_viscosity = @(params) params.mu*1e-6;

mu = params.kinematic_viscosity(params);
c11 = 1e1;% penalty parameter, must be large enough for coercivity
[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 );
required_residual_tol = 1e-6; % allowable residual
max_iter = 2e5; % maximum number of iterations

[ params, paramsP, achieved_residual_tol_schur] =...
    solve_plot_solution_schur( params, paramsP, grid, rhs, stifness_matrix);

[ params, paramsP, flag, achieved_residual_tol, actual_iter] = solve_plot_solution...
    ( params, paramsP, grid, rhs, stifness_matrix, required_residual_tol, max_iter);
 
disp('Entering into stiffness matrix tests')
[ eigen_vectors, eigen_values, condition_number, rank_matrix] = stifness_matrix_test...
( stifness_matrix, params, paramsP, grid, qdeg );

params.dof_analytical = @(glob) [glob(2)*(1-glob(2)) 0];
params.dof_derivative_analytical = @(glob) [0 1-2*glob(y);0 0];
paramsP.dof_analytical = @(glob) (1-glob(1));
paramsP.dof_derivative_analytical = @(glob) [-2*params.kinematic_viscosity(params) 0];
[ error_l2_velocity ] = error_l2_norm_assembly( params, grid, qdeg );
[ error_l2_pressure ] = error_l2_norm_assembly( paramsP, grid, qdeg );

% ERROR FUNCTION CALL

% c11_min = 4e4;
% c11_max = 1e10;
% c11_num_interval = 10;
% [ condition_number, c11 ] = c11_condition_number...
%     ( params, paramsP, grid, qdeg, mu, c11_min, c11_max, c11_num_interval );
% [ solution_norm, c11] = c11_solution( params, paramsP, grid, qdeg,...
%     mu, required_residual_tol, max_iter, c11_min, c11_max, c11_num_interval);