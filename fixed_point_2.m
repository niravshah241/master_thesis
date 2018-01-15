function [ x, iter, update, value ] = fixed_point_2( params, x0, y0, iter_max, tol )
%FIXED_POINT_2 Summary of this function goes here
%   Detailed explanation goes here

%cos(x0)=params.lhs;
%sin(x0)=params.rhs;

for iter=1:1:iter_max
    x = exp(-x0);
    value(iter)=x;
    update(iter)=abs(x-x0);
    if update(iter)<tol
        break
    end
    x0=x;
end
yn=x;
plot(update)
title('update')
figure()
title('value')
plot(value)
pause()
close all
end