function [ res ] = proper_orthogonal_decomposition( params, singular_min, rb_size )
%PROPER_ORTHOGONAL_DECOMPOSITION Summary of this function goes here
%   Detailed explanation goes here

snapshot_matrix = params.snapshots'*params.snapshots;

[U,S,V] = svds(snapshot_matrix,size(snapshot_matrix,1),'largest');
S = diag(real(S));
if params.show_sparsity == true
    figure()
    plot(S)
    title('Singular value decay')
    axis equal
    axis tight
end
rb_size_calculated = sum(singular_min < S);

rb_size = min(rb_size,rb_size_calculated);

rb_basis = U(:,1:rb_size);%check columns or rows are to be taken.
rb_singular_values = S(1:rb_size);

res.rb_basis = rb_basis;
res.rb_singular_values = rb_singular_values;

end