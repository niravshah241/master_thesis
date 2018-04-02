Step_size = [1/4 1/6 1/8 1/10 1/12 1/15 1/17 1/20 1/22 1/25];
Error_velocity = [58.3897 58.4694 55.8762 53.0069 50.3039 46.7447 ...
    44.6924 42.0267 40.4837 38.4248];
Error_pressure = [0.3689 0.3969 0.4111 0.4181 0.4236 0.4293 ...
    0.4320 0.4357 0.4381 0.4410];
figure
plot(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Stokes flow')
axis tight
figure
loglog(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Stokes flow (Logarithmic scale)')
axis tight
figure()
loglog(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Stokes flow (Logarithmic scale)')
axis tight

Step_size = [1/4 1/6 1/8 1/10 1/12 1/15 1/17 1/20 1/22 1/25];
Error_velocity = [326.4299 218.6365 151.7758 110.7830 84.2128 59.1140 ...
    48.0682 36.3153 30.7034 24.3979];
Error_pressure = [0.0186 0.0153 0.0123 0.0102 0.0086 0.0070 ...
    0.0062 0.0053 0.0049 0.0044];
figure()
plot(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Stokes flow')
axis tight
figure
loglog(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Stokes flow (Logarithmic scale)')
axis tight
figure()
loglog(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Stokes flow (Logarithmic scale)')
axis tight