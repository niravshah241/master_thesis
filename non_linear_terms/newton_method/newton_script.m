function [ params,paramsP,flag,relres_solver,iter_solver,...
    relres_newton, iter_newton, stifness_matrix_nonlinear ] =...
    newton_script( params,paramsP,grid,qdeg,mu,c11,...
    tol_newton,max_iter_newton,stifness_matrix, tol_solver, max_iter_solver)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if tol_solver < tol_newton
%     warning('Tolerance for solver should be less than tolerance for newton.')
% end

disp('Entering into solve (newton)');

stifness_matrix_nonlinear = stifness_matrix;

h = zeros(params.ndofs+paramsP.ndofs,1);

for i = 1:1:max_iter_newton
    
    [ res ] = non_linear_term_assembly( params,paramsP,grid,qdeg );
    res_nonlinear = zeros(params.ndofs+paramsP.ndofs);
    res_nonlinear(1:params.ndofs,1:params.ndofs) = res.res;
    stifness_matrix_nonlinear = stifness_matrix_nonlinear + res_nonlinear;
    lhs = 2 * res_nonlinear + stifness_matrix;
    rhs = (- res_nonlinear - stifness_matrix) * [params.dofs;paramsP.dofs] + ...
        [params.linear_side;params.rhs_continuity];
    
    [h, flag, relres_solver, iter_solver] = bicgstab(lhs,rhs,...
        tol_solver,max_iter_solver);
    
    %     [h, flag, relres_solver, iter_solver] = minres(lhs,rhs,...
    %         tol_solver,max_iter_solver);
    
    params.dofs = params.dofs + h(1:params.ndofs,1);
    paramsP.dofs = paramsP.dofs + h(params.ndofs+1:params.ndofs+paramsP.ndofs,1);
    iter_newton = i;
    relres_newton = norm(stifness_matrix_nonlinear * [params.dofs;paramsP.dofs]...
    - [params.linear_side;params.rhs_continuity],2)/norm([params.linear_side;params.rhs_continuity],2);
    
    if relres_newton < tol_newton
        break
    end
    
    close all
end

[ res ] = non_linear_term_assembly( params,paramsP,grid,qdeg );
res_nonlinear = zeros(params.ndofs+paramsP.ndofs);
res_nonlinear(1:params.ndofs,1:params.ndofs) = res.res;

% stifness_matrix_nonlinear = stifness_matrix + res_nonlinear;
stifness_matrix_nonlinear = sparse(stifness_matrix_nonlinear);

relres_newton = norm(stifness_matrix_nonlinear * [params.dofs;paramsP.dofs]...
    - [params.linear_side;params.rhs_continuity],2)/norm([params.linear_side;params.rhs_continuity],2);

disp('entering into plotting Degrees of Freedom (Newton)')

close all

for i=1:1:params.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(params,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Plotting Velocity in x direction (Newton)'])
    else
        title(['Plotting Velocity in y direction (Newton)'])
    end
    %axis equal
    ldg_plot(sdf,grid,params);
    plot(grid);
end

figure()

for i=1:1:paramsP.dimrange
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(paramsP,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure) (Newton)'])
    subplot(paramsP.dimrange,1,i)
    title(['Pressure degree of freedom number ',num2str(i)])
    %axis equal
    ldg_plot(sdf,grid,paramsP);
    plot(grid);
end

end