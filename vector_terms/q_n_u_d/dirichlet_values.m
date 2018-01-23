function [ res ] = dirichlet_values(glob)
%Function evluating dirichlet boundary velocity value based on pointwise evaluation
%   Inputs are global corordinates. Output is column vector

r = norm(glob,2);
x = glob(1);
y = glob(2);
a = 1; %factor for reynolds number


if x>(-eps) & x<(0.1+eps) & y>(1-eps)
    res = a*[10*x 0]';
    
elseif x>0.1 & x<(0.9+eps) & y>(1-eps)
    res = a*[1 0]';
    
elseif x>0.9 & x<=1 & y>(1-eps)
    res = a*[10-10*x 0]';
    
else
    res = a*[0 0]';
    
end

% if glob(1)<eps
%     res = 1*glob(2)*(1-glob(2));
%     res = [res 0]';
%
% elseif (glob(2)>(1-eps) || glob(2)<eps)
%     res = [0 0]';
%
% else
%     res = [0 0]'; %on circular object
% end

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

