function [ res ] = error_h0_norm_assembly( params, grid, qdeg )
%ERROR_L2_NORM_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res = zeros(grid.nelements,1);

for tria_index = 1:1:grid.nelements
    res_integral = error_h0_norm_integral( params, grid, tria_index, qdeg );
    res(tria_index) = res(tria_index) + res_integral^(1/2);
end

res = sum(res);