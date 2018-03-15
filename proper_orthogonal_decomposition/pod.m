function [ pod_res_params, pod_res_paramsP, B_velocity, B_pressure, ...
    red_dim_velocity, red_dim_pressure] = pod( params, paramsP, grid, ...
    red_dim_velocity, red_dim_pressure, min_eigen_pressure, ...
    min_eigen_velocity, stifness_matrix, rhs)
%POD Summary of this function goes here
%   Detailed explanation goes here
% snapshots \in \mathbb(R)^{ndofs \time n_s}

S = params.snapshots_matrix;
n_s = size(S,2);
M_velocity = ldg_mass_matrix(params,grid,params);
snapshots_inner_product = S'*M_velocity*S;

[V,D] = eig(snapshots_inner_product);

[D, iD]=sort(diag(D),'descend');
D=diag(D);
V=V(:,iD);
D = real(D); % remove negligible imaginary part
eigen_values = diag(D);
red_dim_velocity = sum(eigen_values >= min_eigen_velocity);
figure()
plot(eigen_values);
title('Eigen value decay (Velocity) ');
xlabel('Number of eigen value');
ylabel('Eigen values (Velocity)');
R = eye(n_s,red_dim_velocity);
D_sqrt = exp(-0.5*log(eigen_values(1:red_dim_velocity)));
B_velocity = S*V(:,1:red_dim_velocity)*diag(D_sqrt)*R(1:red_dim_velocity,:);
pod_res_params.eigen_values = eigen_values;
pod_res_params.B_velocity = B_velocity;

%Pressure

S = paramsP.snapshots_matrix;
n_s = size(S,2);
M_pressure = ldg_mass_matrix(paramsP,grid,params);
snapshots_inner_product = S'*M_pressure*S;

[V,D] = eig(snapshots_inner_product);

[D, iD]=sort(diag(D),'descend');
D=diag(D);
V=V(:,iD);
D = real(D); % remove negligible imaginary part
eigen_values = diag(D);
red_dim_pressure = sum(eigen_values >= min_eigen_pressure);
figure()
plot(eigen_values);
title('Eigen value decay (Pressure) ');
xlabel('Number of eigen value');
ylabel('Eigen values (Pressure)');
R = eye(n_s,red_dim_pressure);
D_sqrt = exp(-0.5*log(eigen_values(1:red_dim_pressure)));
B_pressure = S*V(:,1:red_dim_pressure)*diag(D_sqrt)*R(1:red_dim_pressure,:);
pod_res_paramsP.eigen_values = eigen_values;
pod_res_paramsP.B = B_pressure;
B_velocity = sparse(B_velocity);
B_pressure = sparse(B_pressure);

end