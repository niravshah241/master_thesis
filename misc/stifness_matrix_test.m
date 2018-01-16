function [ eigen_vectors, eigen_values, condition_number, rank_matrix ] = ...
    stifness_matrix_test( stifness_matrix, params, paramsP, grid, qdeg )
%STIFNESS_MATRIX_TEST Summary of this function goes here
%   Detailed explanation goes here

if ~issparse(stifness_matrix)
    disp('Converting stiffness matrix to sparse matrix (for eigs)')
    stifness_matrix = sparse(stifness_matrix);
end

stifness_matrix = full(stifness_matrix);

A = stifness_matrix - stifness_matrix';

symmetry_zero = find(A~=0);

if (isempty(symmetry_zero))
    disp('Stiffness matrix is symmetric');
else
    disp(['Stiffness matrix is not symmetric with error norm ',...
        num2str(norm(A,2))]);
end

disp('Calculating eigen values and eigen vectors of stiffness matrix, this may take some time')

[V,D] = eigs(stifness_matrix,size(stifness_matrix,1));

D = diag(real(D));

eig_val_non_pos = find(D<=eps);

if (isempty(eig_val_non_pos))
    disp('All eig values are positive');
else
    disp(['There are ',num2str(length(eig_val_non_pos)),' non positive eigen values']);
end

eigen_values = D;
eigen_vectors = V;

disp('Calculating condition number of stiffness matrix, this may take some time')

condition_number = cond(full(stifness_matrix));

disp('Calculating rank of stifness matrix (this may take some time) ')

rank_matrix = rank(full(stifness_matrix));

if nargin>1
    [ res ] = params.lhs_continuity * params.dofs - params.rhs_continuity;
    error = norm(res,2);
    if error == 0
        disp('pressure and velocity fields are orthogonal')
    else
        disp(['pressure and velocity fields are not orthogonal with '...
            'residual ', num2str(error/norm(params.linear_side,2))]);
    end
end
end