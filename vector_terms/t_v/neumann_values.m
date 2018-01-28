function [ res ] = neumann_values(glob,params)
%NEUMANN_VALUES Summary of this function goes here
%   Detailed explanation goes here


% standard condition
if glob(1)>(1-eps)
    res=[0 0]';
else
    disp('Not on neumann boundary');
    %error('Not on neumann boundary');
end


% benchmark problem lid driven cavity
% if glob(2)<eps
%     res=[0 0]';
% else
%     error('Not on neumann boundary');
% end


% analytical example from paper
% if glob(1) < eps
%     res = [0 0]';
% else
%     disp('Not on Neumann boundary')
% end

end