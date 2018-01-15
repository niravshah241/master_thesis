function [ solution_norm, c11] = c11_solution( params, paramsP, grid, qdeg,...
    mu, required_residual_tol, max_iter, c11_min, c11_max, c11_num_interval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

increment = (c11_max-c11_min)/(c11_num_interval);
c11 = c11_min:increment:c11_max;

solution_norm = zeros(length(c11),1);

params_original = params;
paramsP_original = paramsP;

for i=1:1:length(c11)
    
    clear params
    clear paramsP
    clear rhs
    clear stifness_matrix
    
    params = params_original;
    paramsP = paramsP_original;
    
    [ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
        ( params, paramsP, grid, qdeg, mu, c11(i) );
    
    [ params, paramsP, flag, achieved_residual_tol ] = solve_plot_solution...
        ( params, paramsP, grid, rhs, stifness_matrix, required_residual_tol, max_iter);
    
    solution_norm(i) = sqrt(norm(params.dofs,2)^2+norm(paramsP.dofs,2)^2);
    
end

figure()
plot(c11,solution_norm);
title('Solution norm vs c11')

end