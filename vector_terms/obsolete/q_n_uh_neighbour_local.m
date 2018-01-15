function [ res] = q_n_uh_neighbour_local( llcoord, grid, params, tria_index,...
    local_vertex_index)
%Q_N_UH_NEIGHBOUR_LOCAL Summary of this function goes here
%   Detailed explanation goes here
face_index_self = local_vertex_index;
paramsP.pdeg = params.pdeg -1;
paramsP.dimrange=1;
nrep = [3 6 10 15];
paramsP.ndofs_per_element = nrep(paramsP.pdeg)*paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element*grid.nelements;

tria_index_neighbour = grid.NBI(tria_index,face_index_self);
face_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);

pressure_basis_self = ldg_evaluate_basis(llocal2local(grid,...
    face_index_self,llcoord),paramsP);
pressure_basis_neighbour = ldg_evaluate_basis(llocal2local(grid,...
    face_index_neighbour,llcoord),paramsP);
pressure_basis_average = (pressure_basis_neighbour + pressure_basis_self)/2;

velocity_basis_neighbour = ldg_evaluate_basis(llocal2local(grid,...
face_index_neighbour,llcoord),params);

N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
N_neighbour = -N;
jump_velocity_minus = velocity_basis_neighbour*N_neighbour';

for i=1:1:paramsP.ndofs_per_element
    for j=1:1:params.ndofs_per_element    
        res(i,j)=pressure_basis_average(i,:)*jump_velocity_minus(j,:)';
    end
end
end