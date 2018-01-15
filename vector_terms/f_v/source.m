function res=source(lcoord,params, grid, k)

glob=local2global(grid,k,lcoord,params);
f=func_rhs(glob,params);
res=ldg_evaluate_basis(lcoord,params)*f;

end