function [res] = error_velocity_rbasis( params, params_reduced, B_velocity, ...
    grid, qdeg)

res = 0;

for tria_index = 1:1:grid.nelements
    res = res + error_velocity_rbasis_integral(params, params_reduced, ...
    B_velocity, grid, tria_index, qdeg);
end

res = res^0.5;

end