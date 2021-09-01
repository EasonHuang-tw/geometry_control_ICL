classdef controller
    properties
         kx=13;
         kv=10;
         kR = 10*eye(3);
         kW = 1*eye(3);
         %% adaptive
         theta = [0;0;0];
         gamma = 0.0005;
         c2 = 1
        %% ICL
        integral_times_discrete ;
        k_icl;
        N = 10;
        y_i;
        
        y;
        y_omega;
        M_hat;
        
        Y_array;
        Y_omega_array;
        M_array;

        sigma;
         
        end
   methods
              function [control,ex,ev,eR,eW,obj] = geometric_tracking_ctrl(obj,iteration,uav,desired,type)
                desired_p = desired(:,1);
                desired_v = desired(:,2);
                desired_a = desired(:,3);     
                desired_b1= desired(:,4);
                desired_W = [0;0;0];
                desired_W_dot = [0;0;0];
                R_now = reshape(uav.R(:,iteration-1), 3, 3);
                W_now = uav.W(:,iteration-1);
                control = zeros(4,1);
                
                %linear
                ex = uav.x(:,iteration-1) - desired_p;
                ev = uav.v(:,iteration-1) - desired_v;
%                 disp("error_pose")
%                 disp(ex);
                A = (obj.kx*ex + obj.kv*ev + uav.m*uav.g*uav.e3 + desired_a);
                f = dot(A,R_now*uav.e3);
                
                %rotation
                b3_c = A/abs(norm(A));
                b2_c = cross(b3_c,desired_b1);
                b1_c = cross(b2_c,b3_c);
                R_c = [b1_c,b2_c,b3_c];
   
                W_hat = hat(W_now);
                eR = 0.5*vee((R_c'*R_now-R_now'*R_c));     % error R
                eW = W_now-R_now'*R_c*desired_W;

                if type == "origin"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*desired_W - R_now'*R_c*desired_W_dot);
                elseif type == "EMK"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*desired_W - R_now'*R_c*desired_W_dot) - [uav.pc_2_mc(2);-uav.pc_2_mc(1);0]*f;
                elseif type == "adaptive"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*desired_W - R_now'*R_c*desired_W_dot) - obj.theta*f;
                    theta_hat_dot = obj.gamma*f*(eW+obj.c2*eR);
                    obj.theta = obj.theta + theta_hat_dot;
                    obj.theta(3) = 0;
                    disp("theta")
                    disp(obj.theta);
                
                elseif type == "ICL"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*desired_W - R_now'*R_c*desired_W_dot) - obj.theta*f;
                    theta_hat_dot = obj.gamma*f*(eW+obj.c2*eR);
                    obj.theta = obj.theta + theta_hat_dot;
                    obj.theta(3) = 0;
                    
                    obj.y = obj.y + f - obj.Y_array(1);
                    obj.y_omega = obj.y_omega - control.Y_omega_array(:,1);
                    obj.M_hat = obj.M_hat + M -obj.M_array(:,1);
                    
                    for i= 1:obj.integral_times_discrete-1
                        obj.Y_array(i) = obj.Y_array(i+1);
                        obj.Y_omega_array(i) = obj.Y_omega_array(i+1);
                        obj.M_array(i) = obj.M_array(i+1);
                    end
                    if iteration > obj.integral_times_discrete
                    
                    end
                    
                    disp("theta")
                    disp(obj.theta);
                end
%                 disp("M")
%                 disp(M);
                
                control(1) = f;
                control(2:4) = M;
              end
   end
end
