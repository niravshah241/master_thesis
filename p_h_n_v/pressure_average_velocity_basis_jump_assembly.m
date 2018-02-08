function [res] = pressure_average_velocity_basis_jump_assembly(grid,params,paramsP,qdeg)

[ tria_index_internal,local_edge_index_internal,a_index_internal, b_index_internal,...
    local_vertex_index_internal] = tria_edge_index_internal( grid );
[ tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet, b_index_dirichlet,...
    local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

tria_index = [tria_index_internal, tria_index_dirichlet];
local_edge_vertex = [local_edge_index_internal, local_edge_index_dirichlet];
local_vertex_index = [local_vertex_index_internal, local_vertex_index_dirichlet];

gids_velocity = ldg_global_dof_index(params, grid);
gids_pressure = ldg_global_dof_index(paramsP, grid);

% Initialising matrx for each term

res_plus_plus=zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements);
res_plus_minus=zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements);
res_minus_plus=zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements);
res_minus_minus=zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements);

%Assembly of four matrices

for i=1:1:length(tria_index)
    res_plus_plus(gids_pressure(tria_index(i),:),gids_velocity(tria_index(i),:))...
        =res_plus_plus(gids_pressure(tria_index(i),:),gids_velocity(tria_index(i),:))+...
        pressure_average_velocity_basis_jump_int_plus_int_plus...
        (grid,params,paramsP,tria_index(i),local_vertex_index(i),qdeg);
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus(gids_pressure(tria_index(i),:),gids_velocity(tria_index_neighbour,:))...
            =res_plus_minus(gids_pressure(tria_index(i),:),gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_plus_int_minus...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i),qdeg);
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_minus_plus(gids_pressure(tria_index_neighbour,:),gids_velocity(tria_index(i),:))...
            =res_minus_plus(gids_pressure(tria_index_neighbour,:),gids_velocity(tria_index(i),:))+...
            pressure_average_velocity_basis_jump_int_minus_int_plus...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i),qdeg);
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_minus_minus(gids_pressure(tria_index_neighbour,:),gids_velocity(tria_index_neighbour,:))...
            =res_minus_minus(gids_pressure(tria_index_neighbour,:),gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_minus_int_minus...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i),qdeg);
    end
end

% Plotting and storgage in output
close all
res.plus_plus = res_plus_plus;
res.plus_minus = res_plus_minus;
res.minus_plus = res_minus_plus;
res.minus_minus = res_minus_minus;
res.res = res_plus_plus + res_plus_minus + res_minus_plus + res_minus_minus;
if params.show_sparsity == true
    figure()
    spy(res.plus_plus)
    title('Spy of ({\psi^+},[n^+ \cdot \phi^+])')
    figure()
    spy(res.plus_minus)
    title('Spy of ({\psi^+},[n^- \cdot \phi^-])')
    figure()
    spy(res.minus_plus)
    title('Spy of ({\psi^-},[n^+ \cdot \phi^+])')
    figure()
    spy(res.minus_minus)
    title('Spy of ({\psi^-},[n^- \cdot \phi^-])')
    figure()
    spy(res.res)
    title('Spy of ({\psi},[n \cdot \phi])')
    sprintf('Observe all graphs %s','{\psi},[n \cdot \phi]')
    pause();
    close all
end
%3 entries(31,32,33) in diagonal of res.plus_minus and res.minus_plus are non zero
end