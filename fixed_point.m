function [ yn, iter, update, value ] = fixed_point( params, x0, y0, n_max,tol )
%FIXED_POINT Summary of this function goes here
%   Detailed explanation goes here: params should have params.step_size and
%   params.fixed_point_handle fields

%h = params.step_size;
%phi = params.fixed_point_handle;

h=0.1;

phi = @(x,y) fixed_point_handle(x,y);

y_old = y0 + phi(x0,y0)*h;
x = x0+h;

i=1;

for iter=1:1:n_max
    yn = y0 + h * phi(x,y_old);
    value(i) = yn;
    update(i) = abs(yn-y_old);
    if update(i)<tol
        break
    end
    y_old = yn;
    i=i+1;
end
close all
figure()
plot(update)
title(['update'])
figure()
plot(value)
title(['value'])
pause()
close all
end