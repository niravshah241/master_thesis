function [res] = evaluate_jump(llcoord,grid,params,tria_index,local_vertex_index)
%Function evaluating jump of VECTORIAL basis function at internal boundary of grid.
%The basis functions are to be created as row vectors. Output is column vector
face_index = local_vertex_index;
basis = ldg_evaluate_basis(llocal2local(grid,face_index,llcoord),params);
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
face_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
basis_neighbour = ldg_evaluate_basis...
    (llocal2local(grid,face_index_neighbour,llcoord),params);
normal_vector = [grid.NX(tria_index,local_vertex_index)...
    grid.NY(tria_index,local_vertex_index)];
res = (basis - basis_neighbour)*normal_vector';