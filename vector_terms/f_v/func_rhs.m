function [ res ] = func_rhs( x,params )
%FUNC_RHS is function for rhs. res column vector.
%x is global coordinate

res = zeros(params.dimrange,1);

%res = [2*params.kinematic_viscosity(params)-1 0]';

res = [0 0]';
end