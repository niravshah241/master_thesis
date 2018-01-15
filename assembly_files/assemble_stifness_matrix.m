function [ params, paramsP, rhs, stifness_matrix] = assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 )
%variational form assembly
% Momentum : AU + BP = C
% Continuity      B'U = D
params.bilinear_side = assemble_bilinear_side(grid, params,...
    paramsP, qdeg, mu, c11); %A
params.bilinear_side_pressure_terms =...
    assemble_bilinear_side_pressure_terms( params, paramsP, grid, qdeg);%B
params.linear_side = assemble_linear_side(grid, params, qdeg, mu, c11);%C

params.rhs_continuity = assemble_rhs_continuity( params, paramsP, grid, qdeg);%D
params.lhs_continuity = assemble_lhs_continuity( params, paramsP, grid, qdeg);%B'


dofs = [params.dofs;paramsP.dofs];
linear_stifness_matrix = params.linear_side; 
continuity_pressure = zeros(paramsP.ndofs); 
bilinear_stifness_matrix = [params.bilinear_side,params.bilinear_side_pressure_terms];
continuity_stifness_matrix = [params.lhs_continuity,continuity_pressure];
stifness_matrix = [bilinear_stifness_matrix;continuity_stifness_matrix];
rhs = [params.linear_side;params.rhs_continuity];
end

