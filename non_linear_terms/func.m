function [ value, derivative ] = func( x )
%FUNC Summary of this function goes here
%   Detailed explanation goes here

value = x.*sind(x) - cosd(x);
derivative = sind(x) + x.*cosd(x) + sind(x);

end