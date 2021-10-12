clear;
close all;

addpath('geometry-toolbox')
%% set drone parameters
% simulation time
dt = 1/400;
sim_t = 200;

uav = drone_dynamic;
uav.dt = dt;            %delta t
uav.sim_t = sim_t;      %whole 
uav.t = 0:dt:sim_t;     %every time stamps

uav.d = 0.2;            %wing span
uav.m = 1.15;
uav.J = [0.0131, 0, 0;
         0, 0.0131, 0;
         0, 0, 0.0244];

uav.allocation_matrix = cal_allocation_matrix(uav.d, uav.c_tau);
uav.allocation_matrix_inv = cal_allocation_matrix_inv(uav.allocation_matrix);

%% create states array
zero_for_compare = zeros(1,length(uav.t));
uav.x = zeros(3, length(uav.t));
uav.v = zeros(3, length(uav.t));
uav.a = zeros(3, length(uav.t));
uav.x_pose = zeros(3, length(uav.t));
uav.v_pose = zeros(3, length(uav.t));
uav.R = zeros(9, length(uav.t));
uav.W = zeros(3, length(uav.t));
uav.W_dot = zeros(3, length(uav.t));
uav.ex = zeros(3, length(uav.t));
uav.ev = zeros(3, length(uav.t));
uav.eR = zeros(3, length(uav.t));
uav.eW = zeros(3, length(uav.t));
uav.force_moment = zeros(4, length(uav.t));
uav.rotor_thrust = zeros(4, length(uav.t));

real_theta_array = zeros(3, length(uav.t));
theta_array = zeros(3, length(uav.t));

desired_x = zeros(3, length(uav.t));
dX = zeros(18, 1);

%% create Kalmen filter class
kf = z_axis_KF;
kf.r_mp_z = 0;

states_array_x = zeros(3, length(uav.t));
states_array_v = zeros(3, length(uav.t));
states_array_z = zeros(1,length(uav.t));
real_array_z = zeros(1,length(uav.t));
array_y_telta = zeros(6,length(uav.t));
%% initial state
uav.x(:, 1) = [1; 0; 0];
uav.v(:, 1) = [0; 0; 0];
uav.R(:, 1) = [1; 0; 0; 0; 1; 0; 0; 0; 1];
uav.W(:, 1) = [0; 0; 0];

kf.states(1:3,1) = uav.x(:,1);
kf.states(4:6,1) = uav.v(:,1);
kf.states(7,1) =  0;
disp(kf.states);
%% kalmen filter covarience
kf.P =  [   1,    0,    0,    0,    0,    0,    0;
               0, 1,    0,    0,    0,    0,    0;
               0,    0, 1,    0,    0,    0,    0;
               0,    0,    0, 1,    0,    0,    0;
               0,    0,    0,    0, 1,    0,    0;
               0,    0,    0,    0,    0, 1,    0;
               0,    0,    0,    0,    0,    0, 1;];
           
kf.Q = [    0.1,    0,    0,    0,    0,    0,    0;
               0, 0.1,    0,    0,    0,    0,    0;
               0,    0, 0.1,    0,    0,    0,    0;
               0,    0,    0, 0.1,    0,    0,    0;
               0,    0,    0,    0, 0.1,    0,    0;
               0,    0,    0,    0,    0, 0.1,    0;
               0,    0,    0,    0,    0,    0, 0.1;];
           
kf.R = [    0.01,    0,    0,    0,    0,    0;
               0, 0.01,    0,    0,    0,    0;
               0,    0, 0.01,    0,    0,    0;
               0,    0,    0, 0.01,    0,    0;
               0,    0,    0,    0, 0.01,    0;
               0,    0,    0,    0,    0, 0.01;];

%% create controller
control = controller;
integral_time = 0.01;
control.integral_times_discrete = integral_time/uav.dt;
disp(control.integral_times_discrete);
control.y = 0;
control.y_omega = zeros(3,1);
control.M_hat = zeros(3,1);

control.Y_array = zeros(1,control.integral_times_discrete);
control.Y_omega_array = zeros(3,control.integral_times_discrete);
control.M_array = zeros(3,control.integral_times_discrete);
control.W_array = zeros(3,control.integral_times_discrete);

control.sigma_M_hat_array = zeros(3,control.N);
control.sigma_y_omega_array = zeros(3,control.N);
control.sigma_y_array = zeros(control.N);

disp("integral times")
disp(control.integral_times_discrete)
%% create trajectory
traj_array = zeros(12, length(uav.t));
traj = trajectory;

%% create allocation matrix
       allocation_M = cal_allocation_matrix(uav.d,uav.c_tau);
       inv_allocation_M = inv(allocation_M);
       
%% create allocation matrix
        cos_45 = cosd(45);
       uav.pc_2_mc = [0.1;0.1;0.1]; %pose center to mass center
       uav_l = uav.d*cos_45;
       pc_2_r = [  uav_l - uav.pc_2_mc(1),   uav_l - uav.pc_2_mc(1), -(uav_l + uav.pc_2_mc(1)), -(uav_l + uav.pc_2_mc(1));
                   uav_l - uav.pc_2_mc(2),-(uav_l + uav.pc_2_mc(2)), -(uav_l + uav.pc_2_mc(2)),     uav_l- uav.pc_2_mc(2);
                                        0,                        0,                         0,                        0;];
%        disp(pc_2_r(:,1));
%% start iteration

traj_type = "circle";   %"circle","position"
controller_type = "ICL";   %"origin","EMK","adaptive","ICL"

for i = 2:length(uav.t)
    t_now = uav.t(i);
    desired = traj.traj_generate(t_now,traj_type);
    desired_x(:,i) = desired(:,1);
    
    % calculate control force
    [control_output, uav.ex(:, i), uav.ev(:, i), uav.eR(:, i), uav.eW(:, i),control] = control.geometric_tracking_ctrl(i,uav,desired,controller_type);

    % calculate real force applied on the drone
    real_theta_array(:,i) = [uav.pc_2_mc(2),-uav.pc_2_mc(1),0];
    theta_array(:,i) = control.theta;
    rotor_force = allocation_M\ control_output;
    real_control_force = zeros(4,1);
    
    for rotor_num = 1:4
        real_control_force(1) = real_control_force(1)+rotor_force(rotor_num);
        real_control_force(2:4) = real_control_force(2:4) + cross(pc_2_r(:,rotor_num),[0,0,-rotor_force(rotor_num)])';
    end
        real_control_force(4) = [-uav.c_tau, uav.c_tau, -uav.c_tau, uav.c_tau]*rotor_force;

    %% update states
    X0 = [uav.x(:, i-1);
        uav.v(:, i-1);
        reshape(reshape(uav.R(:, i-1), 3, 3), 9, 1);

        uav.W(:, i-1)];
    [T, X_new] = ode45(@(t, x) uav.dynamics( x, real_control_force), [0, dt], X0);

    dX = uav.dynamics(X0 , real_control_force);
    
    % calculate position of pose sensor
    R_wb = reshape(X_new(end, 7:15),3,3);
    x_pose = X_new(end, 1:3) + (R_wb*-uav.pc_2_mc)';
    v_pose = X_new(end, 4:6) + (R_wb*cross(X_new(end, 16:18),-uav.pc_2_mc)')';
    
    %% kalmen filter
    r_pm_x = -control.theta(2,1);
    r_pm_y = control.theta(1,1);
    kf = kf.KF(x_pose,v_pose,r_pm_x,r_pm_y,control_output,uav.m,X_new(end, 1:3),X_new(end, 4:6),X_new(end, 7:15),X_new(end, 16:18),uav.dt);
    disp('kf x');
    disp(kf.states(3));
    disp(X_new(end,3));
    disp(kf.states(3) - X_new(end,3))
    states_array_x(:,i) = kf.states(1:3);
    states_array_v(:,i) = kf.states(4:6);
    states_array_z(1,i) = kf.states(7);
    real_array_z(1,i) = -uav.pc_2_mc(3);
    array_y_telta(:,i) = kf.y_telta;
    
    % Save the states 
    uav.x(:, i) = X_new(end, 1:3);
    uav.v(:, i) = X_new(end, 4:6);
    uav.a(:, i) = dX(4:6)';
    uav.R(:, i) = X_new(end, 7:15);
    uav.W(:, i) = X_new(end, 16:18);
    uav.W_dot(:, i) = dX(16:18)';
    uav.x_pose(:, i) = x_pose;
    uav.v_pose(:, i) = v_pose;
    

end

%% show the result
figure('Name','linear result');

subplot(3,2,1);
plot(uav.t(2:end),uav.x(1,2:end),uav.t(2:end),desired_x(1,2:end));
title('position x')
axis([-inf inf -2 2])
legend({'x','x_d'},'Location','southwest')
subplot(3,2,2);
plot(uav.t(2:end),uav.ex(1,2:end));
title('position error x')
axis([-inf inf -0.5 0.5])

subplot(3,2,3);
plot(uav.t(2:end),uav.x(2,2:end),uav.t(2:end),desired_x(2,2:end));
title('position y')
axis([-inf inf -2 2])
legend({'y','y_d'},'Location','southwest')
subplot(3,2,4);
plot(uav.t(2:end),uav.ex(2,2:end));
title('position error y')
axis([-inf inf -0.5 0.5])

subplot(3,2,5);
plot(uav.t(2:end),uav.x(3,2:end),uav.t(2:end),desired_x(3,2:end));
title('position z')
axis([-inf inf -2 2])
legend({'z','z_d'},'Location','southwest')
subplot(3,2,6);
plot(uav.t(2:end),uav.ex(3,2:end));
title('position error z')
axis([-inf inf -0.5 0.5])
%% rotation
figure('Name','rotation result');

subplot(3,1,1);
plot(uav.t(2:end),uav.eR(1,2:end));
title('eR x')
axis([-inf inf -0.1 0.1])
subplot(3,1,2);
plot(uav.t(2:end),uav.eR(2,2:end));
title('eR y')
axis([-inf inf -0.1 0.1])
subplot(3,1,3);
plot(uav.t(2:end),uav.eR(3,2:end));
title('eR z')
axis([-inf inf -0.1 0.1])

%% theta
figure('Name','theta result');

subplot(2,1,1);
plot(uav.t(2:end), theta_array(1,2:end),uav.t(2:end),real_theta_array(1,2:end));
legend({'theta1','theta1_d'},'Location','southwest')
title('theta 1')
axis([-inf inf 0.05 0.15])
subplot(2,1,2);
plot(uav.t(2:end), theta_array(2,2:end),uav.t(2:end),real_theta_array(2,2:end));
legend({'theta2','theta2_d'},'Location','southwest')
title('theta 2')
axis([-inf inf -0.15 -0.05])

%% theta
figure('Name','pose');

subplot(3,1,1);
plot(uav.t(2:end), uav.x(1,2:end),uav.t(2:end),uav.x_pose(1,2:end));
legend({'x','x_pose'},'Location','southwest')
title('x')
axis([-inf inf -2 2])
subplot(3,1,2);
plot(uav.t(2:end), uav.x(2,2:end),uav.t(2:end),uav.x_pose(2,2:end));
legend({'y','y_pose'},'Location','southwest')
title('y')
axis([-inf inf -2 2])
subplot(3,1,3);
plot(uav.t(2:end), uav.x(3,2:end),uav.t(2:end),uav.x_pose(3,2:end));
legend({'z','z_pose'},'Location','southwest')
title('z')
axis([-inf inf -2 2])

%% position and velocity of position sensor
figure('Name','vel');

subplot(3,1,1);
plot(uav.t(2:end), uav.v(1,2:end),uav.t(2:end),uav.v_pose(1,2:end));
legend({'x_vel','x_p_vel'},'Location','southwest')
title('x')
axis([-inf inf -2 2])
subplot(3,1,2);
plot(uav.t(2:end), uav.v(2,2:end),uav.t(2:end),uav.v_pose(2,2:end));
legend({'y_vel','y_p_vel'},'Location','southwest')
title('y')
axis([-inf inf -2 2])
subplot(3,1,3);
plot(uav.t(2:end), uav.v(3,2:end),uav.t(2:end),uav.v_pose(3,2:end));
legend({'z_vel','z_p_vel'},'Location','southwest')
title('z')
axis([-inf inf -2 2])

%% position and velocity of position sensor
% figure('Name','omega');
% 
% subplot(3,1,1);
% plot(uav.t(2:end), uav.W(1,2:end));
% title('x')
% axis([-inf inf -2 2])
% subplot(3,1,2);
% plot(uav.t(2:end), uav.W(2,2:end));
% title('y')
% axis([-inf inf -2 2])
% subplot(3,1,3);
% plot(uav.t(2:end), uav.W(3,2:end));
% title('z')
% axis([-inf inf -2 2])

%% states
figure('Name','KF_result_x');

subplot(3,1,1);
plot(uav.t(2:end), states_array_x(1,2:end),uav.t(2:end),uav.x(1,2:end));
legend({'x_measured','x_real'},'Location','southwest')
title('x')
axis([-inf inf -2 2])
subplot(3,1,2);
plot(uav.t(2:end), states_array_x(2,2:end),uav.t(2:end),uav.x(2,2:end));
legend({'y_measured','y_real'},'Location','southwest')
title('y')
axis([-inf inf -2 2])
subplot(3,1,3);
plot(uav.t(2:end), states_array_x(3,2:end),uav.t(2:end),uav.x(3,2:end));
legend({'z_measured','z_real'},'Location','southwest')
title('z')
axis([-inf inf -2 2])

figure('Name','KF_result_v');
subplot(3,1,1);
plot(uav.t(2:end), states_array_v(1,2:end),uav.t(2:end),uav.v(1,2:end));
legend({'xv_measured','xv_real'},'Location','southwest')
title('xv')
axis([-inf inf -4 4])
subplot(3,1,2);
plot(uav.t(2:end), states_array_v(2,2:end),uav.t(2:end),uav.v(2,2:end));
legend({'yv_measured','yv_real'},'Location','southwest')
title('yv')
axis([-inf inf -4 4])
subplot(3,1,3);
plot(uav.t(2:end), states_array_v(3,2:end),uav.t(2:end),uav.v(3,2:end));
legend({'zv_measured','zv_real'},'Location','southwest')
title('zv')
axis([-inf inf -4 4])


disp(length( states_array_z));
disp(length(uav.t));
figure('Name','KF_result_z');
plot(uav.t(2:end), states_array_z(1,2:end),uav.t(2:end), real_array_z(1,2:end));
legend({'zv_measured','zv_real'},'Location','southwest')
title('zv')
axis([-inf inf -1 1])

%% states error
figure('Name','KF_result_error');

subplot(4,1,1);
plot(uav.t(2:end), states_array_x(1,2:end)-uav.x(1,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'x_error'},'Location','southwest')
title('x')
axis([-inf inf -0.5 0.5])
subplot(4,1,2);
plot(uav.t(2:end), states_array_x(2,2:end)-uav.x(2,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'y_error'},'Location','southwest')
title('y')
axis([-inf inf -0.5 0.5])
subplot(4,1,3);
plot(uav.t(2:end), states_array_x(3,2:end)-uav.x(3,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'z_error'},'Location','southwest')
title('z')
axis([-inf inf -0.5 0.5])
subplot(4,1,4);
plot(uav.t(2:end), states_array_z(1,2:end) - real_array_z(1,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'r_pm_z_error'},'Location','southwest')
title('r_pm_z')
axis([-inf inf -0.5 0.5])

figure('Name','KF_vel_result_error');

subplot(3,1,1);
plot(uav.t(2:end), states_array_v(1,2:end)-uav.v(1,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'x_error'},'Location','southwest')
title('x')
axis([-inf inf -0.5 0.5])
subplot(3,1,2);
plot(uav.t(2:end), states_array_v(2,2:end)-uav.v(2,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'y_error'},'Location','southwest')
title('y')
axis([-inf inf -0.5 0.5])
subplot(3,1,3);
plot(uav.t(2:end), states_array_v(3,2:end)-uav.v(3,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'z_error'},'Location','southwest')
title('z')
axis([-inf inf -0.5 0.5])

%% disp y_telta
figure('Name','KF_result_error');

subplot(6,1,1);
plot(uav.t(2:end), array_y_telta(1,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'x'},'Location','southwest')
title('x')
axis([-inf inf -0.001 0.001])
subplot(6,1,2);
plot(uav.t(2:end), array_y_telta(2,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'y_error'},'Location','southwest')
title('y')
axis([-inf inf -0.001 0.001])
subplot(6,1,3);
plot(uav.t(2:end), array_y_telta(3,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'z_error'},'Location','southwest')
title('z')
axis([-inf inf -0.001 0.001])
subplot(6,1,4);
plot(uav.t(2:end), array_y_telta(4,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'vx'},'Location','southwest')
title('vx')
axis([-inf inf -0.001 0.001])
subplot(6,1,5);
plot(uav.t(2:end), array_y_telta(5,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'vy'},'Location','southwest')
title('vy')
axis([-inf inf -0.001 0.001])
subplot(6,1,6);
plot(uav.t(2:end), array_y_telta(6,2:end),uav.t(2:end),zero_for_compare(1,2:end));
legend({'vz'},'Location','southwest')
title('vz')
axis([-inf inf -0.001 0.001])