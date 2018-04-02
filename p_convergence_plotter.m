figure()
error_l2_velocity = load('elvelocity.mat');
error_l2_velocity = error_l2_velocity.error_l2_velocity ;
h = load('pconvstepsize.mat');
h = h.h;
pde = [2 3 4];
plot(h,error_l2_velocity(:,1),'-.o',h,error_l2_velocity(:,2),'-.o',h,error_l2_velocity(:,3),'-.o')
xlabel('Step size')
ylabel('Velocity L^2 error')
title('p-convergence for velocity for Stokes equation (L^2 error)')
legend('D = 2','D = 3', 'D = 4')
axis tight

figure()
error_l2_pressure = load('elpressure.mat');
error_l2_pressure = error_l2_pressure.error_l2_pressure ;
step_size = h;
pde = [2 3 4];
plot(h,error_l2_pressure(:,1),'-.o',h,error_l2_pressure(:,2),'-.o',h,error_l2_pressure(:,3),'-.o')
xlabel('Step size')
ylabel('Pressure L^2 error')
title('p-convergence for pressure for Stokes equation (L^2 error)')
legend('D-1 = 1','D-1 = 2', 'D-1 = 3')
axis tight

% figure()
% error_h0_pressure = load('ehpressure.mat');
% error_h0_pressure = error_h0_pressure.error_h0_pressure ;
% step_size = h;
% pde = [2 3 4];
% plot(h,error_h0_pressure(:,1),h,error_h0_pressure(:,2),h,error_h0_pressure(:,3))
% xlabel('Step size')
% ylabel('Pressure H_0 error')
% title('p-convergence for pressure for Stokes equation (H_0 error)')
% legend('D = 2','D = 3', 'D = 4')
% axis tight

figure()
error_h0_velocity = load('ehvelocity.mat');
error_h0_velocity = error_h0_velocity.error_h0_velocity ;
step_size = h;
pde = [2 3 4];
plot(h,error_h0_velocity(:,1),'-.o',h,error_h0_velocity(:,2),'-.o',h,error_h0_velocity(:,3),'-.o')
xlabel('Step size')
ylabel('Velocity H_0 error')
title('p-convergence for velocity for Stokes equation (H_0 error)')
legend('D = 2','D = 3', 'D = 4')
axis tight