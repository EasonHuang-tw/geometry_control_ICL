classdef drone_dynamic
    properties
        % simulation time
        dt = 0.01;
        sim_t = 5;     %total simulation duration
        t;              %ever time steps
        iter
        % parameters
        m = 1.25
        J = eye(3)
        d = 0.2         %wing span
        pc_2_mc
        c_tau = 1.347e-2;
        g = 9.81;
        allocation_matrix
        allocation_matrix_inv
        % unit vector
        e1 = [1; 0; 0];
        e2 = [0; 1; 0];
        e3 = [0; 0; 1];
        % states
        x
        v
        a
        R
        W
        W_dot
        
        x_pose
        v_pose
        % errors
        ex
        ev
        eR
        eW
        % control input
        force_moment = zeros(1,4)
        rotor_thrust = zeros(1,4)
    end
    
    methods
        function dX = dynamics(obj , X, F)
            dX = zeros(18, 1);  %delta of pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,Rotation_matrix*9, w_x,w_y,w_z
            R_now = reshape(X(7:15), 3, 3);
            W_now = X(16:18);
            f = F(1);
            M = F(2:4);            
            dx = X(4:6);
            dv = obj.g*obj.e3 - (f/obj.m)*R_now*obj.e3;
            dR = R_now*hat_map(W_now);
            dW = obj.J\(-cross(W_now, obj.J*W_now) + M);
            
            dX(1:3) = dx;
            dX(4:6) = dv;
            dX(7:15) = reshape(dR, 9, 1);
            dX(16:18) = dW;
        end
    end
end