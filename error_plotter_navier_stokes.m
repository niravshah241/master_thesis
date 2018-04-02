Step_size = [1/6 1/8 1/10 1/12 1/15 1/17 1/20];
Error_velocity = [63.2463, 57.9821, 49.8227, 44.7326, 41.4134, 39.7150, ...
    37.2217];
Error_pressure = [3.6204, 4.4499, 3.7357, 3.5708, 4.1540, 4.6418, ...
    5.2330];
figure
plot(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Navier Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Navier Stokes flow')
axis tight
figure
loglog(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('H_0 Error Velocity')
title('Step size vs H_0 Error in Velocity for Navier Stokes flow (Logarithmic scale)')
axis tight
figure()
loglog(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('H_0 Error Pressure')
title('Step size vs H_0 Error in Pressure for Navier Stokes flow (Logarithmic scale)')
axis tight

Step_size = [1/6 1/8 1/10 1/12 1/15 1/17 1/20];
Error_velocity = [173.8074, 127.2059, 76.8123, 56.0014, 40.1646, 32.3621, ...
    25.0632];
Error_pressure = [0.7543, 0.8544, 0.4964, 0.3904, 0.4395, 0.4937, ...
    0.5531];
figure()
plot(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Navier Stokes flow')
axis tight
figure()
plot(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Navier Stokes flow')
axis tight
figure
loglog(Step_size,Error_velocity,'-.or')
xlabel('Step size')
ylabel('L^2 Error Velocity')
title('Step size vs L^2 Error in Velocity for Navier Stokes flow (Logarithmic scale)')
axis tight
figure()
loglog(Step_size,Error_pressure,'-.or')
xlabel('Step size')
ylabel('L^2 Error Pressure')
title('Step size vs L^2 Error in Pressure for Navier Stokes flow (Logarithmic scale)')
axis tight