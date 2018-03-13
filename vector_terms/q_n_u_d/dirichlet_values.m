function [ res ] = dirichlet_values(glob,params)
%Function evluating dirichlet boundary velocity value based on pointwise evaluation
%   Inputs are global corordinates. Output is column vector

if ~isfield(params,'parameter_training_set')
    
    r = norm(glob,2);
    x = glob(1);
    y = glob(2);
    
    %
    %% analytical example from paper
    % if y < eps || y > (1-eps) || x > (1-eps)
    %     res = [0 0]';
    % else
    %     res = [0 0]';
    % end
    
    %% benchmark problem
    % a = 1e2; %factor for reynolds number
    % if x>(-eps) & x<(0.1+eps) & y>(1-eps)
    %     res = a*[10*x 0]';
    %
    % elseif x>0.1 & x<(0.9+eps) & y>(1-eps)
    %     res = a*[1 0]';
    %
    % elseif x>0.9 & x<=1 & y>(1-eps)
    %     res = a*[10-10*x 0]';
    %
    % else
    %     res = a*[0 0]';
    %
    % end
    
    %% standard
    if glob(1)<eps
        res = 1e0*glob(2)*(1-glob(2));
        res = [res 0]';
        
    elseif (glob(2)>(1-eps) || glob(2)<eps)
        res = [0 0]';
        
    else
        res = [0 0]'; %on circular object
    end
else
    r = norm(glob,2);
    x = glob(1);
    y = glob(2);
    
    %
    %% analytical example from paper
    % if y < eps || y > (1-eps) || x > (1-eps)
    %     res = params.parameter_training_set(2)*[0 0]';
    % else
    %     res = params.parameter_training_set(2)*[0 0]';
    % end
    
    %% benchmark problem
    % a = params.parameter_training_set(2); %factor for reynolds number
    % if x>(-eps) & x<(0.1+eps) & y>(1-eps)
    %     res = a*[10*x 0]';
    %
    % elseif x>0.1 & x<(0.9+eps) & y>(1-eps)
    %     res = a*[1 0]';
    %
    % elseif x>0.9 & x<=1 & y>(1-eps)
    %     res = a*[10-10*x 0]';
    %
    % else
    %     res = a*[0 0]';
    %
    % end
    
    
    %% standard
    if glob(1)<eps
        res = -params.parameter_training_set(2)*1e0*glob(2)*(1-glob(2));
        res = [res 0]';
        
    elseif (glob(2)>(1-eps) || glob(2)<eps)
        res = params.parameter_training_set(2)*[0 0]';
        
    else
        res = params.parameter_training_set(2)*[0 0]'; %on circular object
    end
end

end