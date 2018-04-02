function [params, paramsP] = dg_solution( params, paramsP, grid, ...
    reference_factor, qdeg)


params.parameter_training_set = params.parameter_online;

tic;

[ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, params.parameter_online(1), ...
    reference_factor * params.parameter_online(1));

time_matrix_assembly_affine = toc;
disp(['Time taken for assembling stifness matrix ',...
    num2str(time_matrix_assembly_affine)])


[ params, paramsP, achieved_residual_tol_schur] = solve_plot_solution_schur...
    ( params, paramsP, grid, rhs, stifness_matrix);

for i=1:1:params.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(params,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (Schur)'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Velocity in x direction (Schur)'])
    else
        title(['Velocity in y direction (Schur)'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,grid,params);
    plot(grid);
end

for i=1:1:paramsP.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(paramsP,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure)'])
    %subplot(paramsP.dimrange,1,i)
    title('Pressure (Schur)')
    axis equal
    axis tight
    ldg_plot(sdf,grid,paramsP);
    plot(grid);
end
end