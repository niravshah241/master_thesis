function [ res ] = assemble_rhs_continuity( params, paramsP, grid, qdeg )
%ASSEMBLE_ Summary of this function goes here
%   Detailed explanation goes here

rhs = q_n_u_d_assembly( params, paramsP, grid, qdeg );
%rhs = zeros(size(rhs));
rhs = sparse(rhs);

res=rhs;

end

