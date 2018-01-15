function [ res ] = c_h_u_v_term_1_integral( params, paramsP, grid, tria_index, qdeg )
%C_H_U_V_TERM_1_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(lcoord) c_h_u_v_term_1_local( lcoord, params, paramsP, grid, tria_index);

res = -triaquadrature(qdeg,f)*grid.A(tria_index);
%grid.A(tria_index) = detrminant of jacobian

end

