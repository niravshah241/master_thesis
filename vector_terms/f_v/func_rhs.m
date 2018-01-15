function [ res ] = func_rhs( x,params )
%FUNC_RHS is function for rhs. res column vector.
%x is global coordinate
res=zeros(params.dimrange,1);
res(1:length(x))=...
    [(-3*(x(1))^2/(params.kinematic_viscosity(params)))+...
    (2*params.kinematic_viscosity(params))*x(2) 0];

end

