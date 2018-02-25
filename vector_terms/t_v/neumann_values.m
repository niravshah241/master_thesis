function [ res ] = neumann_values(glob,params)
%NEUMANN_VALUES Summary of this function goes here
%   Detailed explanation goes here

%% analytical example from paper
% if glob(1) < eps
%     res = [0 0]';
% else
%     disp('Not on Neumann boundary')
% end

%% benchmark problem lid driven cavity
b = params.parameter_training_set(3);
if glob(2)<eps
    res=b*[0 0]';
else
    error('Not on neumann boundary');
end

%% standard condition
% if glob(1)>(1-eps)
%     res=[0 0]';
% else
%     disp('Not on neumann boundary');
%     %error('Not on neumann boundary');
% end

end