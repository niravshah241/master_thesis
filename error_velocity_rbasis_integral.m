function [res] = error_velocity_rbasis_integral(params, params_reduced, ...
    B_velocity, grid, tria_index, qdeg)

f = @(lcoord) error_velocity_rbasis_local(lcoord, params, params_reduced, ...
    B_velocity, grid, tria_index);

res = triaquadrature(qdeg,f) * 2 * grid.A(tria_index);