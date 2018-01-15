function res=source_integral(params, grid, k, qdeg);
%B=[grid.JIT(k,:,1);grid.JIT(k,:,2)];
f=@(lcoord) source(lcoord, params, grid, k);
res=triaquadrature(qdeg,f)*2*grid.A(k);%2*grid.A is jacobian determinant
end