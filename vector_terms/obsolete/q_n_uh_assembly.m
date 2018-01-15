function [res] = q_n_uh_assembly(grid, params, qdeg);

paramsP.pdeg = params.pdeg -1;
paramsP.dimrange=1;
nrep = [3 6 10 15];
paramsP.ndofs_per_element = nrep(paramsP.pdeg)*paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element*grid.nelements;

[ tria_index,local_edge_index,a_index, b_index, local_vertex_index] =...
    tria_edge_index_internal( grid );

res_plus = zeros(paramsP.ndofs,params.ndofs);
res_minus = zeros(paramsP.ndofs,params.ndofs);

gids_velocity = ldg_global_dof_index(params,grid);
gids_pressure = ldg_global_dof_index(paramsP,grid);


for i=1:1:length(tria_index)
    res_plus(gids_pressure(tria_index(i),:),gids_velocity(tria_index(i),:))...
        =res_plus(gids_pressure(tria_index(i),:),gids_velocity(tria_index(i),:))...
        + q_n_uh_self_integral( grid, params, tria_index(i),local_vertex_index(i), qdeg);
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    res_minus(gids_pressure(tria_index(i),:),gids_velocity(tria_index_neighbour,:))...
        =res_minus(gids_pressure(tria_index(i),:),gids_velocity(tria_index_neighbour,:))...
        + q_n_uh_neighbour_integral( grid, params, tria_index(i),local_vertex_index(i), qdeg);
end

res.res = res_plus + res_minus;
res.res_plus = res_plus;
res.res_minus = res_minus;

figure()
spy(res.res_plus)
figure()
spy(res.res_minus)
%TODO: some res.res_minus entries are non-zero
end