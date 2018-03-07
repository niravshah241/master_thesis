function [ pod_res_params, pod_res_paramsP, B, B_velocity, B_pressure, ...
    stifness_matrix_reduced, params_reduced, paramsP_reduced] = ...
    pod( params, paramsP, grid, red_dim_velocity, red_dim_pressure, ...
    min_eigen, stifness_matrix, rhs, error_norm)
%POD Summary of this function goes here
%   Detailed explanation goes here
% snapshots \in \mathbb(R)^{ndofs \time n_s}

S = params.snapshots_matrix;
n_s = size(S,2);
M = ldg_mass_matrix(params,grid,params);
M_velocity = M;
snapshots_inner_product = S'*M*S;

[V,D] = eig(snapshots_inner_product);

[D, iD]=sort(diag(D),'descend');
D=diag(D);
V=V(:,iD);

D = real(D); % remove negligible imaginary part

eigen_values = diag(D);
figure()
plot(eigen_values);
title('Eigen value decay (Velocity) ');
xlabel('Number of eigen value');
ylabel('Eigen values (Velocity)');
R = eye(n_s,red_dim_velocity);
B_velocity = S*V*D^(-0.5)*R;
pod_res_params.eigen_values = eigen_values;
pod_res_params.B_velocity = B_velocity;

%Pressure

S = paramsP.snapshots_matrix;
n_s = size(S,2);
M = ldg_mass_matrix(paramsP,grid,params);
M_pressure = M;
snapshots_inner_product = S'*M*S;

[V,D] = eig(snapshots_inner_product);

[D, iD]=sort(diag(D),'descend');
D=diag(D);
V=V(:,iD);

D = real(D); % remove negligible imaginary part

eigen_values = diag(D);
figure()
plot(eigen_values);
title('Eigen value decay (Pressure) ');
xlabel('Number of eigen value');
ylabel('Eigen values (Pressure)');
R = eye(n_s,red_dim_pressure);
B_pressure = S*V*D^(-0.5)*R;
pod_res_paramsP.eigen_values = eigen_values;
pod_res_paramsP.B = B_pressure;

% Galerkin formulation
B = zeros(params.ndofs+paramsP.ndofs,red_dim_velocity+red_dim_pressure);
B(1:params.ndofs,1:red_dim_velocity) = B_velocity;
B(params.ndofs+1:params.ndofs+paramsP.ndofs,...
    red_dim_velocity+1:red_dim_velocity+red_dim_pressure) = B_pressure;
B(red_dim_velocity+1:red_dim_velocity+red_dim_pressure,...
    params.ndofs+1:params.ndofs+paramsP.ndofs) = B_pressure';
stifness_matrix_reduced = B'*stifness_matrix*B;
rhs_reduced = B' * rhs;
if error_norm == 'L2'
    params_reduced.dofs = B_velocity' * M_velocity * params.dofs;
    paramsP_reduced.dofs = B_pressure' * M_pressure * paramsP.dofs;
elseif error_norm == 'energy'
    reduced_dofs = minres(stifness_matrix_reduced,rhs_reduced);
    params_reduced.dofs = reduced_dofs(1:red_dim_velocity,1);
    paramsP_reduced.dofs = reduced_dofs...
        (red_dim_velocity+1:red_dim_velocity+red_dim_pressure,1);
    params_reduced.dofs = reduced_dofs(1:red_dim_velocity,1);
end
end