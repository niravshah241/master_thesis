Step_size = [1/10 1/15 1/17 1/20];
Error_velocity = [46.8707 38.6485 36.2001 33.1710];
Error_pressure = [2.2698 0.3407 0.3369 0.3320];
figure
plot(Step_size,Error_velocity)
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Navier Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure)
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Navier Stokes flow')
axis tight
figure
plot(log(Step_size),log(Error_velocity))
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Navier Stokes flow (Logarithmic scale)')
axis tight
figure()
plot(log(Step_size),log(Error_pressure))
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Navier Stokes flow (Logarithmic scale)')
axis tight

Step_size = [1/10 1/15 1/17 1/20];
Error_velocity = [75.7447 36.9755 28.2105 19.8837];
Error_pressure = [0.2436 0.0067 0.0056 0.0051];
figure
plot(Step_size,Error_velocity)
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Navier Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure)
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Navier Stokes flow')
axis tight
figure
plot(log(Step_size),log(Error_velocity))
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Navier Stokes flow (Logarithmic scale)')
axis tight
figure()
plot(log(Step_size),log(Error_pressure))
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Navier Stokes flow (Logarithmic scale)')
axis tight