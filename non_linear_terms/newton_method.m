function [ value_newton, x, update, value, derivative, iter ] =...
    newton_method( x0, nmax, tol_residual, tol)
%NEWTON_METHOD Summary of this function goes here
%   Detailed explanation goes here

[ value(1), derivative(1) ] = func( x0 );

x(1) = x0;

for i = 1:1:nmax
    x(i+1) = x(i) - value(i)./derivative(i);
    update(i) = norm(x(i+1)-x(i),2);
    if  update(i)< tol
        disp(['The newton method stopped due to update tolerance criteria ',...
            'the last update of solution is ',newline, num2str(update(i)),...
            ' is less update tolerance ',num2str(tol)])
        break
    end
    [ value(i+1), derivative(i+1) ] = func( x(i+1) );
    if abs(value(i+1)) < tol_residual
        disp(['The newton method stopped due to residual tolerance criteria ',...
            'the last residual of solution is ',newline, num2str(value(i+1)),...
            ' is less than residual tolerance ',num2str(tol_residual)])
        break
    end
    if i == nmax
        disp(['The number of iterations ',num2str(i),...
            ' reached maximum number of allowable iteration ',newline, num2str(nmax),...
            ' The residual tolerance = ', num2str(value(i+1))])
        
    end
end

iter = i;
value_newton = value(length(value));

end