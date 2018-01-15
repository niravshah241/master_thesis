function [res] = pressure_average_velocity_basis_jump_local_minus_local_plus...
    (llcoord, grid,params,paramsP,tria_index,local_vertex_index)
%function evaluating L^2 norm
%(\psi^-, n^+ \cdot \phi^+ )_{\partial \Omega \cup \partial \Omega_D}
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
face_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
N = [grid.NX(tria_index, local_vertex_index)...
    grid.NY(tria_index, local_vertex_index)];
% N is normal
q = ldg_evaluate_basis(llocal2local(grid,face_index_neighbour,llcoord),paramsP);
% q is pressure basis
v = ldg_evaluate_basis(llocal2local(grid,local_vertex_index,llcoord),params);
%v is velocity basis and face_index_self = local_vertex_index
v_dot_n = v*N';
res = zeros(paramsP.ndofs_per_element,params.ndofs_per_element);

if tria_index_neighbour>0
    
        for i=1:1:paramsP.ndofs_per_element
            for j=1:1:params.ndofs_per_element
                res(i,j) = res(i,j)+q(i)*v_dot_n(j);
            end
        end
    
else
    %do nothing as \psi^- does not exist on dirichlet boundary
end
end