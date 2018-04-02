function [params, paramsP, params_reduced, paramsP_reduced] = online_phase...
    (params, paramsP, grid, stifness_matrix_reference, rhs_reference, ...
    B_velocity, B_pressure, linear_side_reference, red_dim_velocity, ...
    red_dim_pressure, reference_factor)

stifness_matrix = stifness_matrix_reference;
rhs = rhs_reference;
theta_1 = params.parameter_online(1) / params.reference_parameter(1);
%viscocity
theta_2 = params.parameter_online(2) / params.reference_parameter(2);
%dirichlet value
stifness_matrix = stifness_matrix_reference;
rhs = rhs_reference;
disp('Assembling stifness matrix affine')
tic
stifness_matrix(1:params.ndofs,1:params.ndofs) = ...
    theta_1 * stifness_matrix(1:params.ndofs,1:params.ndofs);

rhs(1:params.ndofs) = rhs(1:params.ndofs) + ...
    (params.parameter_online(1) * theta_2 - ...
    params.reference_parameter(1)) * (reference_factor * ...
    linear_side_reference.term_3 - linear_side_reference.term_4);

rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
    theta_2 * rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);

time_matrix_assembly_affine = toc;
disp(['Time taken for affine assembling stifness matrix ',...
    num2str(time_matrix_assembly_affine)])

[params_reduced, paramsP_reduced, stifness_matrix_reduced, rhs_reduced] = ...
    galerkin_formulation(stifness_matrix, rhs, params, paramsP, grid, ...
    red_dim_velocity, red_dim_pressure, B_velocity, B_pressure);

velocity_reduced = B_velocity * params_reduced.dofs;
pressure_reduced = B_pressure * paramsP_reduced.dofs;
params_reduced.dofs = velocity_reduced;
paramsP_reduced.dofs = pressure_reduced;

params_reduced.dimrange = params.dimrange;
params_reduced.ndofs_per_element = params.ndofs_per_element;
params_reduced.ndofs = params.ndofs;
paramsP_reduced.dimrange = paramsP.dimrange;
paramsP_reduced.ndofs_per_element = paramsP.ndofs_per_element;
paramsP_reduced.ndofs = paramsP.ndofs;

for i=1:1:params_reduced.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(params_reduced,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (Schur Reduced)'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Velocity in x direction (Schur Reduced)'])
    else
        title(['Velocity in y direction (Schur Reduced)'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,grid,params);
    plot(grid);
end

for i=1:1:paramsP_reduced.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(paramsP_reduced,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure Reduced)'])
    %subplot(paramsP.dimrange,1,i)
    title('Pressure (Schur Reduced)')
    axis equal
    axis tight
    ldg_plot(sdf,grid,paramsP);
    plot(grid);
end

% error_velocity.dofs = params.dofs - params_reduced.dofs;
% error_velocity.pdeg = params.pdeg;
% error_velocity.dimrange = params.dimrange;
% error_velocity.ndofs_per_element = params.ndofs_per_element;
% error_velocity.ndofs = params.ndofs;
% error_velocity.grid = grid;
% error_pressure.dimrange = paramsP.dimrange;
% error_pressure.ndofs_per_element = paramsP.ndofs_per_element;
% error_pressure.ndofs = paramsP.ndofs;
% error_pressure.dofs = paramsP.dofs - paramsP_reduced.dofs;
% error_pressure.pdeg = paramsP.pdeg;
% error_pressure.grid = grid;
% for i=1:1:error_velocity.dimrange
%     figure()
%     axis equal
%     [scalar_dofs, scalar_df_info] = ldg_scalar_component(error_velocity,i);
%     sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
%     disp(['Plotting',num2str(i),' degree of freedom (error velocity)'])
%     %subplot(params.dimrange,1,i)
%     %title(['Velocity degree of freedom number ',num2str(i)])
%     if i==1
%         title(['Error velocity in x direction'])
%     else
%         title(['Error velocity in y direction'])
%     end
%     axis equal
%     axis tight
%     ldg_plot(sdf,grid,error_velocity);
%     plot(grid);
% end
% 
% for i=1:1:error_pressure.dimrange
%     figure()
%     [scalar_dofs, scalar_df_info] = ldg_scalar_component(error_pressure,i);
%     sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
%     disp(['Plotting ',num2str(i),' degree of freedom (Error pressure)'])
%     %subplot(paramsP.dimrange,1,i)
%     title('Error in pressure')
%     axis equal
%     axis tight
%     ldg_plot(sdf,grid,error_pressure);
%     plot(grid);
% end


end