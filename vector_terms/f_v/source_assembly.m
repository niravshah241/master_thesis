function [ res ] = source_assembly(params, grid,qdeg )
%SOURCE_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs,1);
gids = ldg_global_dof_index(params, grid);

for k=1:1:grid.nelements
    A = source_integral(params, grid, k, qdeg);
    res(gids(k,:)) = res(gids(k,:)) + A;
end
end