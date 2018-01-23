function [ res ] = neumann_values(glob,params)
%NEUMANN_VALUES Summary of this function goes here
%   Detailed explanation goes here

if glob(1)>(1-eps)
    res=[0 0]';
else
    disp('Not on neumann boundary');
    %error('Not on neumann boundary');
end

% if glob(2)<eps
%     res=[0 0]';
% else
%     error('Not on neumann boundary');
% end


end