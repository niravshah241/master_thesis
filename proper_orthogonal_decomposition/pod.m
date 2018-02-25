function [ res ] = pod( params, grid, red_dim, min_eigen,...
    stifness_matrix)
%POD Summary of this function goes here
%   Detailed explanation goes here
% snapshots \in \mathbb(R)^{ndofs \time n_s}

S = params.snapshots_matrix;
n_s = size(S,2);
M = ldg_mass_matrix(params,grid,params);

snapshots_inner_product = S'*M*S;

[V,D] = eig(snapshots_inner_product);

[D, iD]=sort(diag(D),'descend'); 
D=diag(D);                      
V=V(:,iD);                      

D = real(D); % remove negligible imaginary part

eigen_values = diag(D);
figure()
plot(eigen_values);
title('Eigen value decay ');
xlabel('Number of eigen value');
ylabel('Eigen values');
R = eye(n_s,red_dim);
B = S*V*D^(-0.5)*R;
% Galerkin formulation
res.eigen_values = eigen_values;
res.B = B;

end