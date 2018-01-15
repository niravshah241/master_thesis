function [ res ] = dirichlet_values(glob)
%Function evluating dirichlet boundary velocity value based on pointwise evaluation
%   Inputs are global corordinates. Output is column vector

r = norm(glob,2);

if glob(1)<eps
    res=1e0*glob(2)*(1-glob(2));
    res=[res 0]';
    
elseif (glob(2)>(1-eps) || glob(2)<eps)
    res=[0 0]';
    
else
    res = [0 0]'; %on circular object   
end

% if glob(2)>(1-eps)
%     res=-100*glob(1)*(1-glob(1));
%     res=[res 0]';
%     
% elseif (glob(1)>(1-eps) || glob(1)<eps)
%     res=[0 0]';
%     
% else
%     res = [0 0]'; %on circular object
% end


end

