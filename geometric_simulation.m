clear;
close all;

addpath('geometry-toolbox')
%% set drone parameters
% simulation time
dt = 1/400;
sim_t =30;

uav1 = drone_dynamic;
uav1.dt = dt;            %delta t
uav1.sim_t = sim_t;      %whole 
uav1.t = 0:dt:sim_t;     %every time stamps

uav1.d = 0.2;            %wing span
uav1.m = 1.15;
uav1.J = [0.05, -0.005, -0.006;
         -0.005, 0.06, 0.004;
         -0.006, 0.004, 0.07];
uav1.allocation_matrix = cal_allocation_matrix(uav1.d, uav1.c_tau);
uav1.allocation_matrix_inv = cal_allocation_matrix_inv(uav1.allocation_matrix);

%% create states array
zero_for_compare = zeros(1,length(uav1.t));
uav1.x = zeros(3, length(uav1.t));
uav1.v = zeros(3, length(uav1.t));
uav1.a = zeros(3, length(uav1.t));
uav1.x_pose = zeros(3, length(uav1.t));
uav1.v_pose = zeros(3, length(uav1.t));
uav1.R = zeros(9, length(uav1.t));
uav1.W = zeros(3, length(uav1.t));
uav1.W_dot = zeros(3, length(uav1.t));
uav1.ex = zeros(3, length(uav1.t));
uav1.ev = zeros(3, length(uav1.t));
uav1.eR = zeros(3, length(uav1.t));
uav1.eW = zeros(3, length(uav1.t));
uav1.force_moment = zeros(4, length(uav1.t));
uav1.rotor_thrust = zeros(4, length(uav1.t));


real_theta_array_uav1 = zeros(8, length(uav1.t));
theta_array_uav1 = zeros(8, length(uav1.t));
theta_hat_dot_array_uav1 = zeros(8, length(uav1.t));
contorl_output_array_uav1 = zeros(4, length(uav1.t));
dX_uav1 = zeros(18, 1);

theta_array_uav2 = zeros(8, length(uav1.t)); %use adaptive

desired_x = zeros(3, length(uav1.t));

%% initial state
uav1.x(:, 1) = [1; 0; 0];
uav1.v(:, 1) = [0; 0; 0];
uav1.R(:, 1) = [1; 0; 0; 0; 1; 0; 0; 0; 1];
uav1.W(:, 1) = [0; 0; 0];


%% create controller
control_uav1 = controller;


control_uav1.sigma_M_hat_array = zeros(3,control_uav1.N);
control_uav1.sigma_y_omega_array = zeros(3,control_uav1.N);
control_uav1.sigma_y_array = zeros(3,8,control_uav1.N);

control_uav2 = control_uav1;
%% create trajectory
traj_array = zeros(12, length(uav1.t));
traj = trajectory;

%% create allocation matrix
       allocation_M = cal_allocation_matrix(uav1.d,uav1.c_tau);
       inv_allocation_M = inv(allocation_M);
       
%% create CoG allocation matrix
        cos_45 = cosd(45);
       uav1.pc_2_mc = [0.1;0.1;0.0]; %pose center to mass center
%        uav.pc_2_mc = [0;0;0];
       uav_l = uav1.d*cos_45;
       mc_2_r = [  uav_l - uav1.pc_2_mc(1),   uav_l - uav1.pc_2_mc(1), -(uav_l + uav1.pc_2_mc(1)), -(uav_l + uav1.pc_2_mc(1));
                   uav_l - uav1.pc_2_mc(2),-(uav_l + uav1.pc_2_mc(2)), -(uav_l + uav1.pc_2_mc(2)),     uav_l- uav1.pc_2_mc(2);
                                        0,                        0,                         0,                        0;];

uav2 = uav1;
              
 %% start iteration

traj_type = "circle";   %"circle","position"
controller_type = "ICL";   %"origin","EMK","adaptive","ICL"

control_output_uav1  = zeros(4,1);
control_output_uav2  = zeros(4,1);
% control_uav1.gamma =  diag([0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00001,0.00001])*0.5;
control_uav2.gamma =  diag([0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00001,0.00001])*0.5;
for i = 2:length(uav1.t)
    disp(i)
    t_now = uav1.t(i);
    desired = traj.traj_generate(t_now,traj_type);
    desired_x(:,i) = desired(:,1);
    
    % calculate control force
    [control_output_uav1, uav1.ex(:, i), uav1.ev(:, i), uav1.eR(:, i), uav1.eW(:, i),control_uav1] = control_uav1.geometric_tracking_ctrl(i,uav1,desired,controller_type);
    [control_output_uav2, uav2.ex(:, i), uav2.ev(:, i), uav2.eR(:, i), uav2.eW(:, i),control_uav2] = control_uav2.geometric_tracking_ctrl(i,uav2,desired,'adaptive');

    % calculate real force applied on the drone
    real_theta_array_uav1(:,i) = [uav1.J(1,1),uav1.J(2,2),uav1.J(3,3),uav1.J(1,2),uav1.J(1,3),uav1.J(2,3),uav1.pc_2_mc(1),uav1.pc_2_mc(2)];
    theta_hat_dot_array_uav1(:,i) =control_uav1.theta_hat_dot;
    
    theta_array_uav1(:,i) = control_uav1.theta;
    theta_array_uav2(:,i) = control_uav2.theta;

    
    rotor_force_uav1 = allocation_M\ control_output_uav1;
    real_control_force_uav1 = zeros(4,1);
    contorl_output_array_uav1(:,i) = control_output_uav1;
    for rotor_num = 1:4
        real_control_force_uav1(1) = real_control_force_uav1(1)+rotor_force_uav1(rotor_num);
        real_control_force_uav1(2:4) = real_control_force_uav1(2:4) + cross(mc_2_r(:,rotor_num),[0,0,-rotor_force_uav1(rotor_num)])';
    end
        real_control_force_uav1(4) = [-uav1.c_tau, uav1.c_tau, -uav1.c_tau, uav1.c_tau]*rotor_force_uav1;
   
    rotor_force_uav2 = allocation_M\ control_output_uav2;
    real_control_force_uav2 = zeros(4,1);
    for rotor_num = 1:4
        real_control_force_uav2(1) = real_control_force_uav2(1)+rotor_force_uav2(rotor_num);
        real_control_force_uav2(2:4) = real_control_force_uav2(2:4) + cross(mc_2_r(:,rotor_num),[0,0,-rotor_force_uav2(rotor_num)])';
    end
        real_control_force_uav2(4) = [-uav2.c_tau, uav2.c_tau, -uav2.c_tau, uav2.c_tau]*rotor_force_uav2;
    %% update states
    X0_uav1 = [uav1.x(:, i-1);
        uav1.v(:, i-1);
        reshape(reshape(uav1.R(:, i-1), 3, 3), 9, 1);
        uav1.W(:, i-1)];
   

    [T_uav1, X_new_uav1] = ode45(@(t, x) uav1.dynamics( x, real_control_force_uav1), [0, dt], X0_uav1);
    
    M_tot_uav1 = real_control_force_uav1(2:4);
    dX_uav1 = uav1.dynamics(X0_uav1 , real_control_force_uav1);
    
     
    % calculate position of pose sensor
    R_wb_uav1 = reshape(X_new_uav1(end, 7:15),3,3);
    x_pose_uav1 = X_new_uav1(end, 1:3) + (R_wb_uav1*-uav1.pc_2_mc)';
    v_pose_uav1 = X_new_uav1(end, 4:6) + (R_wb_uav1*cross(X_new_uav1(end, 16:18),-uav1.pc_2_mc)')';
   
    % Save the states 
    uav1.x(:, i) = X_new_uav1(end, 1:3);
    uav1.v(:, i) = X_new_uav1(end, 4:6);
    uav1.a(:, i) = dX_uav1(4:6)';
    uav1.R(:, i) = X_new_uav1(end, 7:15);
    uav1.W(:, i) = X_new_uav1(end, 16:18);
    uav1.W_dot(:, i) = dX_uav1(16:18)';
    uav1.x_pose(:, i) = x_pose_uav1;
    uav1.v_pose(:, i) = v_pose_uav1;
  %% uav2 update
    X0_uav2 = [uav2.x(:, i-1);
        uav2.v(:, i-1);
        reshape(reshape(uav2.R(:, i-1), 3, 3), 9, 1);
        uav2.W(:, i-1)];
    
    
    [T_uav2, X_new_uav2] = ode45(@(t, x) uav2.dynamics( x, real_control_force_uav2), [0, dt], X0_uav2); 
    
    M_tot_uav2 = real_control_force_uav2(2:4);
    dX_uav2 = uav1.dynamics(X0_uav2 , real_control_force_uav2);
    
    R_wb_uav2 = reshape(X_new_uav2(end, 7:15),3,3);
    x_pose_uav2 = X_new_uav2(end, 1:3) + (R_wb_uav2*-uav2.pc_2_mc)';
    v_pose_uav2 = X_new_uav2(end, 4:6) + (R_wb_uav2*cross(X_new_uav2(end, 16:18),-uav2.pc_2_mc)')';
   
    % Save the states 
    uav2.x(:, i) = X_new_uav2(end, 1:3);
    uav2.v(:, i) = X_new_uav2(end, 4:6);
    uav2.a(:, i) = dX_uav2(4:6)';
    uav2.R(:, i) = X_new_uav2(end, 7:15);
    uav2.W(:, i) = X_new_uav2(end, 16:18);
    uav2.W_dot(:, i) = dX_uav2(16:18)';
    uav2.x_pose(:, i) = x_pose_uav2;
    uav2.v_pose(:, i) = v_pose_uav2;
end

%% show the result
figure('Name','linear result');

subplot(3,2,1);
plot(uav1.t(2:end),uav1.x(1,2:end),"LineWidth",2);
hold on;
plot(uav1.t(2:end),desired_x(1,2:end),'--',"LineWidth",2);
title('$Position\ in\ X\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
ylabel(' $X[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -2 2])
legend('$X$','$X_d$','interpreter','latex','FontSize',24);
subplot(3,2,2);
plot(uav1.t(2:end),uav1.ex(1,2:end),"LineWidth",2);
title('$Position\ Error\ in\ X\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
ylabel(' $eX[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -0.5 0.5])

subplot(3,2,3);
plot(uav1.t(2:end),uav1.x(2,2:end),"LineWidth",2);
hold on;
plot(uav1.t(2:end),desired_x(2,2:end),'--',"LineWidth",2);
title('$Position\ in\ Y\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
ylabel(' $Y[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -2 2])
legend('$Y$','$Y_d$','interpreter','latex','FontSize',24);
subplot(3,2,4);
plot(uav1.t(2:end),uav1.ex(2,2:end),"LineWidth",2);
title('$Position\ Error\ in\ Y\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
ylabel(' $eY[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -0.5 0.5])

subplot(3,2,5);
plot(uav1.t(2:end),uav1.x(3,2:end),"LineWidth",2);
hold on;
plot(uav1.t(2:end),desired_x(3,2:end),'--',"LineWidth",2);
title('$Position\ in\ Z\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
ylabel(' $Z[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -2 2])
legend('$Z$','$Z_d$','interpreter','latex','FontSize',24);
subplot(3,2,6);
plot(uav1.t(2:end),uav1.ex(3,2:end),"LineWidth",2);
title('$Position\ Error\ in\ Z\ Direction\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
ylabel(' $eZ[m]$','interpreter','latex','FontSize',24);
axis([-inf inf -0.5 0.5])
%% rotation
figure('Name','rotation result');

subplot(3,1,1);
plot(uav1.t(2:end),uav1.eR(1,2:end),"LineWidth",2);
legend('$row\ error$','interpreter','latex','FontSize',24);
ylabel(' $eR_x[degree]$','interpreter','latex','FontSize',24);
title('$Attitude\ Errors\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
axis([-inf inf -1 1])
subplot(3,1,2);
plot(uav1.t(2:end),uav1.eR(2,2:end),"LineWidth",2);
legend('$pitch\ error$','interpreter','latex','FontSize',24);
ylabel(' $eR_y[degree]$','interpreter','latex','FontSize',24);
axis([-inf inf -1 1])
subplot(3,1,3);
plot(uav1.t(2:end),uav1.eR(3,2:end),"LineWidth",2);
legend('$yaw\ error$','interpreter','latex','FontSize',24);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
ylabel(' $eR_z[degree]$','interpreter','latex','FontSize',24);
axis([-inf inf -1 1])

%% theta inertial
figure('Name','inertia estimation result');

subplot(6,1,1);
plot(uav1.t(2:end), theta_array_uav1(1,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(1,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,xx}$','$\theta_{J,xx}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,xx}[kgm^2]$','interpreter','latex','FontSize',12);
title('$Estimation\ of\ the\ Moment\ of\ Inertia\ in\ Simulation$','interpreter','latex', 'FontSize', 24);
axis([-inf inf real_theta_array_uav1(1,2)-0.05 real_theta_array_uav1(1,2)+0.05])

subplot(6,1,2);
plot(uav1.t(2:end), theta_array_uav1(2,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(2,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,yy}$','$\theta_{J,yy}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,yy}[kgm^2]$','interpreter','latex','FontSize',12);
axis([-inf inf real_theta_array_uav1(2,2)-0.05 real_theta_array_uav1(2,2)+0.05])

subplot(6,1,3);
plot(uav1.t(2:end), theta_array_uav1(3,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(3,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,zz}$','$\theta_{J,zz}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,zz}[kgm^2]$','interpreter','latex','FontSize',12);
axis([-inf inf real_theta_array_uav1(3,2)-0.05 real_theta_array_uav1(3,2)+0.05])

subplot(6,1,4);
plot(uav1.t(2:end), theta_array_uav1(4,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(4,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,xy}$','$\theta_{J,xy}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,xy}[kgm^2]$','interpreter','latex','FontSize',12);
axis([-inf inf real_theta_array_uav1(4,2)-0.01 real_theta_array_uav1(4,2)+0.01])

subplot(6,1,5);
plot(uav1.t(2:end), theta_array_uav1(5,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(5,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,xz}$','$\theta_{J,xz}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,xz}[kgm^2]$','interpreter','latex','FontSize',12);
axis([-inf inf real_theta_array_uav1(5,2)-0.01 real_theta_array_uav1(5,2)+0.01])

subplot(6,1,6);
plot(uav1.t(2:end), theta_array_uav1(6,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(6,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{J,yz}$','$\theta_{J,yz}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{J,yz}[kgm^2]$','interpreter','latex','FontSize',12);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
axis([-inf inf real_theta_array_uav1(6,2)-0.01 real_theta_array_uav1(6,2)+0.01])

%% theta CoG
figure('Name','CoG estimation result');
subplot(2,1,1);
plot(uav1.t(2:end), theta_array_uav1(7,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(7,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{{r_{mp}}_x}$','$\theta_{{r_{mp}}_x}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{{r_{mp}}_x}[m]$','interpreter','latex','FontSize',24);
title('$Estimation\ of\ the\ Center\ of\ Gravity\ in\ Simulation$','interpreter','latex', 'FontSize', 24);

subplot(2,1,2);
plot(uav1.t(2:end), theta_array_uav1(8,2:end),"LineWidth",2,'Color',[0 0.5 1]);
hold on;
plot(uav1.t(2:end),real_theta_array_uav1(8,2:end),'--',"LineWidth",2,'Color',[1 0.3 0]);
hold off;
legend('$\hat{\theta}_{{r_{mp}}_y}$','$\theta_{{r_{mp}}_y}$','interpreter','latex','FontSize',24);
ylabel(' $\hat{\theta}_{{r_{mp}}_y}[m]$','interpreter','latex','FontSize',24);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
%% theta_hat_dot
% figure('Name','theta_hat_dot result');
% 
% subplot(8,1,1);
% plot(uav1.t(2:end), theta_hat_dot_array_uav1(1,2:end));
% legend({'theta hat dot 1'},'Location','southwest')
% title('theta hat dot 1')
% 
% subplot(8,1,2);
% plot(uav1.t(2:end), theta_hat_dot_array_uav1(2,2:end));
% legend({'theta hat dot 2'},'Location','southwest')
% title('theta hat dot 2')
% 
% subplot(8,1,3);
% plot(uav1.t(2:end), theta_hat_dot_array_uav1(3,2:end));
% legend({'theta hat dot 3'},'Location','southwest')
% title('theta hat dot 3')
% 
% subplot(8,1,4);
% plot(uav1.t(2:end), theta_hat_dot_array_uav1(4,2:end));
% legend({'theta hat dot 4'},'Location','southwest')
% title('theta hat dot 4')

%% theta_norm_compare
compare_array_uav1 = zeros(1,length(uav1.t));
for i = 1:length(uav1.t)
    compare_array_uav1(i) = norm(theta_array_uav1(:,i)-real_theta_array_uav1(:,2))/norm(real_theta_array_uav1(:,2));
end

compare_array_uav2 = zeros(1,length(uav1.t));
for i = 1:length(uav1.t)
    compare_array_uav2(i) = norm(theta_array_uav2(:,i)-real_theta_array_uav1(:,2))/norm(real_theta_array_uav1(:,2));
end

figure('Name','theta_norm_compare');
plot(uav1.t(2:end),compare_array_uav1(2:end),"LineWidth",2);
hold on;
plot(uav1.t(2:end),compare_array_uav2(2:end),"LineWidth",2);
hold off;
legend('$\frac{\|\tilde{\theta}\|_{adaptive}}{\|\theta\|}$','$\frac{\|\tilde{\theta}\|_{ICL}}{\|\theta\|}$','interpreter','latex','FontSize',24);
set(gca, 'YScale', 'log');
ylabel(' $\frac{\|\tilde{\theta}\|}{\|\theta\|}$','interpreter','latex','FontSize',24);
xlabel(' $t[s]$','interpreter','latex','FontSize',24);
title('$Estimation\ Errors\ of\ Adaptive\ and\ ICL$','interpreter','latex', 'FontSize', 24);

%% contorl output
% figure('Name','control output result');
% 
% subplot(4,1,1);
% plot(uav1.t(2:end), contorl_output_array_uav1(1,2:end));
% legend({'control 1'},'Location','southwest')
% title('control 1')
% 
% subplot(4,1,2);
% plot(uav1.t(2:end), contorl_output_array_uav1(2,2:end));
% legend({'control 2'},'Location','southwest')
% title('control 2')
% 
% subplot(4,1,3);
% plot(uav1.t(2:end), contorl_output_array_uav1(3,2:end));
% legend({'control 3'},'Location','southwest')
% title('control 3')
% 
% subplot(4,1,4);
% plot(uav1.t(2:end), contorl_output_array_uav1(3,2:end));
% legend({'control 4'},'Location','southwest')
% title('control 4')
