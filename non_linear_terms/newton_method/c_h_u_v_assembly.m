function [ res ] = c_h_u_v_assembly( params, paramsP, grid, qdeg )
%C_H_U_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);

res_volume_integral = zeros(params.ndofs);

for tria_index = 1:1:grid.nelements
    temp = c_h_u_v_term_1_integral( params, paramsP, grid, tria_index, qdeg );
    res_volume_integral(gids(tria_index,:),gids(tria_index,:)) = ...
        res_volume_integral(gids(tria_index,:),gids(tria_index,:)) + temp;
end
res.volume_integral = 2*res_volume_integral;

if params.show_sparsity == true
    figure()
    spy(res.volume_integral)
    title('Spy of newton method 2(-w \cdot \phi,\phi)')
end

res_term_4 = zeros(params.ndofs);

[tria_index_neumann,local_edge_index_neumann,a_index_neumann,...
    b_index_neumann, local_vertex_index_neumann] = tria_edge_index_neumann( grid );

for i = 1:1:length(tria_index_neumann)
    temp = c_h_u_v_term_4_integral( params,...
        paramsP, grid, tria_index_neumann(i), local_vertex_index_neumann(i), qdeg);
    res_term_4(gids(tria_index_neumann(i),:),gids(tria_index_neumann(i),:)) = ...
        res_term_4(gids(tria_index_neumann(i),:),gids(tria_index_neumann(i),:)) + temp;
end
res.term_4 = 2*res_term_4;

if params.show_sparsity == true
    figure()
    spy(res.term_4)
    title('Spy of newton method 2((w \cdot n) \phi^T \phi)')
end

[tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet,...
    b_index_dirichlet, local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

[tria_index_internal,local_edge_index_internal,a_index_internal,...
    b_index_internal, local_vertex_index_internal] = tria_edge_index_internal( grid );

tria_index = [tria_index_dirichlet,tria_index_internal];
local_vertex_index = [local_vertex_index_dirichlet,local_vertex_index_internal];

res_term_2_1 = zeros(params.ndofs);

for i=1:1:length(tria_index)
    temp = c_h_u_v_term_2_1_integral( params, paramsP,...
        grid, tria_index(i), local_vertex_index(i), qdeg);
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_term_2_1(gids(tria_index(i),:),gids(tria_index_neighbour,:)) = ...
            res_term_2_1(gids(tria_index(i),:),gids(tria_index_neighbour,:)) + temp;
    end
end
res.term_2_1 = 2*res_term_2_1;

if params.show_sparsity == true
    figure()
    spy(res.term_2_1)
    title('Spy of newton method 2*1/2*((w \cdot n_i) \phi {\phi^{ext}}^T)')
end

res_term_2_2 = zeros(params.ndofs);

for i=1:1:length(tria_index)
    temp = c_h_u_v_term_2_2_integral( params, paramsP,...
        grid, tria_index(i), local_vertex_index(i), qdeg);
    res_term_2_2(gids(tria_index(i),:),gids(tria_index(i),:)) = ...
        res_term_2_2(gids(tria_index(i),:),gids(tria_index(i),:)) + temp;
end
res.term_2_2 = 2*res_term_2_2;

if params.show_sparsity == true
    figure()
    spy(res.term_2_2)
    title('Spy of newton method 2*1/2*((w \cdot n_i) \phi {\phi}^T)')
end

res_term_3_1 = zeros(params.ndofs);

for i=1:1:length(tria_index)
    temp = c_h_u_v_term_3_1_integral( params, paramsP,...
        grid, tria_index(i), local_vertex_index(i), qdeg);
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_term_3_1(gids(tria_index(i),:),gids(tria_index_neighbour,:)) = ...
            res_term_3_1(gids(tria_index(i),:),gids(tria_index_neighbour,:)) + temp;
    end
end
res.term_3_1 = 2*res_term_3_1;

if params.show_sparsity == true
    figure()
    spy(res.term_3_1)
    title('Spy of newton method 2*1/2*(abs(w \cdot n_i) \phi {\phi^{ext}}^T)')
end

res_term_3_2 = zeros(params.ndofs);

for i=1:1:length(tria_index)
    temp = c_h_u_v_term_2_2_integral( params, paramsP,...
        grid, tria_index(i), local_vertex_index(i), qdeg);
    res_term_3_2(gids(tria_index(i),:),gids(tria_index(i),:)) = ...
        res_term_3_2(gids(tria_index(i),:),gids(tria_index(i),:)) + temp;
end
res.term_3_2 = 2*res_term_3_2;

if params.show_sparsity == true
    figure()
    spy(res.term_3_2)
    title('Spy of newton method 2*1/2*(abs(w \cdot n_i) \phi {\phi}^T)')
end

res.res = res.volume_integral + res.term_2_1 + res.term_2_2 - res.term_3_1 +...
    res.term_3_2 + res.term_4;

end