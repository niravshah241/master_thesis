function [ res ] = non_linear_term_assembly( params,paramsP,grid,qdeg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res.res1 = zeros(params.ndofs);
res.res2 = zeros(params.ndofs);
res.res3 = zeros(params.ndofs);
res.res4 = zeros(params.ndofs);
res.res = zeros(params.ndofs);

[ res1 ] = w_dot_del_u_v_assembly(params, paramsP, grid, qdeg );
[ res2 ] = w_ni_u_v_assembly( params, paramsP, grid, qdeg );
[ res3 ] = abs_w_ni_u_v_assembly( params, paramsP, grid, qdeg );
[ res4 ] = w_n_u_v_assembly( params, paramsP, grid, qdeg);

res.res1 = sparse(res1);
res.res2 = sparse(res2.res);
res.res3 = sparse(res3.res);
res.res4 = sparse(res4);
res.res = sparse(res1 + res2.res + res3.res + res4);

if params.show_sparsity == true
    figure()
    spy(full(res.res))
    title('spy of non linear term c(u;u,\phi (self part))')
end

end