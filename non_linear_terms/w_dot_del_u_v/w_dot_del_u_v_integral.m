function [ res] = w_dot_del_u_v_integral( params, paramsP, grid, tria_index, qdeg)
%W_DOT_DEL_U_V_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(lcoord) w_dot_del_u_v_local( lcoord, params, paramsP, grid, tria_index);

res = -triaquadrature(qdeg,f)*2*grid.A(tria_index); 

end