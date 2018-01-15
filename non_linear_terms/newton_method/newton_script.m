function [ params,paramsP,flag,relres,iter ] = newton_script( params,paramsP,grid,qdeg,mu,c11, tol_solution,max_non_linear_iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = zeros(params.ndofs,1);
% TODO 1: NEEDS FULL CORRECTION in formulation
for i=1:1:max_non_linear_iter
    [ res1 ] = c_h_u_v_assembly( params, paramsP, grid, qdeg );%2*a
    [ res2 ] = non_linear_term_assembly( params,paramsP,grid,qdeg );%a
    [ res3 ] = assemble_linear_side( grid,params,qdeg, mu, c11);%f
    
    right_hand_side = sparse(res2.res*params.dofs-res3);
    left_hand_side = sparse(res1.res);
    
    [h,flag,relres,iter] = bicgstab(left_hand_side,right_hand_side,tol_solution,max_non_linear_iter);
    %TODO 2: suitable solver and more input output arguments i.e. diferent
    %tol_solution and max_non_linear_iter than for for loop
    
    params.dofs = params.dofs + h;
    
    if norm(h,2)<tol_solution
        break
    end
end    
end