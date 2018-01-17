function [ params,paramsP,flag,relres_solver,iter_solver,...
    relres_newton, iter_newton, stifness_matrix_nonlinear ] =...
    newton_script( params,paramsP,grid,qdeg,mu,c11,...
    tol_newton,max_iter_newton,stifness_matrix, tol_solver, max_iter_solver)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

stifness_matrix_nonlinear = stifness_matrix;

h = zeros(params.ndofs+paramsP.ndofs,1);

for i = 1:1:max_iter_newton
    
    [ res ] = non_linear_term_assembly( params,paramsP,grid,qdeg );
    res_nonlinear = zeros(params.ndofs+paramsP.ndofs);
    res_nonlinear(1:params.ndofs,1:params.ndofs) = res;
    stifness_matrix_nonlinear = stifness_matrix_nonlinear + res_nonlinear;
    lhs = 2 * res_nonlinear + stifness_matrix;
    rhs = (- res_nonlinear - stifness_matrix) * [params.dofs;paramsP.dofs] + ...
        [params.linear_side;params.rhs_continuity];
    
    %1. check formula
    
    [h, flag, relres_solver, iter_solver] = bicgstab(lhs,rhs,...
        tol_solver,max_iter_solver);
    
    % [h, flag, relres_solver, iter_solver] = minres(lhs,rhs,...
    %     tol_solver,max_iter_solver);
    
    params.dofs = params.dofs + h(1:params.ndofs,1);
    paramsP.dofs = paramsP.dofs + h(params.ndofs+1:params.ndofs+paramsP.ndofs,1);
    iter_newton = i;
    
    if norm(h,2) < tol_newton
        break
    end
      
end

[ res ] = non_linear_term_assembly( params,paramsP,grid,qdeg );
res_nonlinear = zeros(params.ndofs+paramsP.ndofs);
res_nonlinear(1:params.ndofs,1:params.ndofs) = res;
stifness_matrix_nonlinear = stifness_matrix_nonlinear + res_nonlinear;
relres_newton = stifness_matrix_nonlinear * [params.dofs;paramsP.dofs]...
    - [params.linear_side;params.rhs_continuity];
% check relres formula
% check sizes of differenet matrices
end