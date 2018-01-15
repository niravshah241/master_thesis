function [ res ] = newton_gradient( params, gids, k )
%NEWTON_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.dimrange);

res = [cos(params.dofs(gids(k,:))) sin(params.dofs(gids(k,:)))];

end