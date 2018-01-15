function [ res ] = newton_func( params,gids,k)
%NEWTON_FUNC Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.dimrange,1);

res = sin(params.dofs(gids(k,:))) - cos(params.dofs(gids(k,:)));

end

