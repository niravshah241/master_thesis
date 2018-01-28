function [ res ] = error_h0_norm_integral( params, grid, tria_index, qdeg )
%ERROR_H0_NORM_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(lcoord) error_h0_norm_local( lcoord, params, grid, tria_index);

res = triaquadrature(qdeg,f)*2*grid.A(tria_index);

end