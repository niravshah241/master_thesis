function [ res ] = fixed_point_handle( x,y )
%FIXED_POINT_HANDLE Summary of this function goes here
%   Detailed explanation goes here

for i=1:1:size(x,1)
    %input x is such that x(i,:) represents one entry i.e. x is arranged row wise)
    
    res(i) = exp(sin(x(i,:)))*(cos(x(i,:)));
    
end

end

