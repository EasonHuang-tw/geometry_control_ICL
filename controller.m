classdef controller
    properties
         kx=13;
         kv=6;
         kR = 14*eye(3);
         kW = 1.2*eye(3);
         %% adaptive
         theta = [0;0;0];
         gamma =  0.00007;
%         gamma = 0;
         c2 = 12
        %% ICL
        integral_times_discrete ;
        k_icl = 0.005;
        N = 10;
        y_i;
        
        y;
        y_omega;
        M_hat;
        
        Y_array;
        Y_omega_array;
        M_array;
        W_array;
        R_array;

        sigma_y_array;
        sigma_y_omega_array;
        sigma_M_hat_array;
         
        end
   methods
             function [u, u_dot, u_ddot] = deriv_unit_vector(obj,q, q_dot, q_ddot)

                nq = norm(q);
                u = - q / nq;
                u_dot = - q_dot / nq + dot(q, q_dot)*q / nq^3;

                u_ddot = - q_ddot / nq +  (2 * dot(q, q_dot))*q_dot / nq^3 ...
                    +  (dot(q_dot, q_dot) + dot(q, q_ddot))*q / nq^3 ...
                    - 3 * dot(q, q_dot)^2 * q/ nq^5;

              end
              function [control,ex,ev,eR,eW,obj] = geometric_tracking_ctrl(obj,iteration,uav,desired,type)
                desired_p = desired(:,1);
                desired_v = desired(:,2);
                desired_a = desired(:,3);     
                desired_j = desired(:,4);  
                desired_s = desired(:,5);
                desired_b1= desired(:,6);
                desired_b1_dot= desired(:,7);
                desired_b1_ddot= desired(:,8);
                desired_W = [0;0;0];
                desired_W_dot = [0;0;0];
                R_now = reshape(uav.R(:,iteration-1), 3, 3);
                W_now = uav.W(:,iteration-1);
                control = zeros(4,1);
                
                %linear
                ex = uav.x(:,iteration-1) - desired_p;
                ev = uav.v(:,iteration-1) - desired_v;
                ea = uav.a(:,iteration-1) - desired_a;
%                 disp("error_pose")
%                 disp(ex);
                

                A = -obj.kx*ex - obj.kv*ev - uav.m*uav.g*uav.e3 + uav.m*desired_a;
                f = dot(-A,R_now*uav.e3);
                %% calc A_dot,A_ddot
                
                
                A_dot = -obj.kx * ev - obj.kv * ea + uav.m * desired_j;
                
                b3 = R_now * uav.e3;
                b3_dot = R_now * hat(W_now) * uav.e3;
                f_dot = -dot(A_dot, b3) - dot(A, b3_dot);
                ej = - f_dot / uav.m * b3 - f / uav.m * b3_dot - desired_j;
                
                A_ddot = -obj.kx * ea - obj.kv * ej + uav.m * desired_s;
                
                [b3_c, b3_c_dot, b3_c_ddot] = obj.deriv_unit_vector(A, A_dot, A_ddot);
                %disp(A_dot);

                A2 = -hat(b3_c) * desired_b1;
                A2_dot = -(hat(b3_c_dot) * desired_b1 + hat(b3_c) * desired_b1_dot);
                A2_ddot = -(hat(b3_c_ddot) * desired_b1 + 2 * hat(b3_c_dot) * desired_b1_dot + hat(b3_c) * desired_b1_ddot);
                [b2_c, b2_c_dot, b2_c_ddot] = obj.deriv_unit_vector(A2, A2_dot, A2_ddot);

                %rotation
                %b3_c = -A/abs(norm(A));
                %b2_c = cross(b3_c,desired_b1);
                b1_c = hat(b2_c) * b3_c;
                b1_c_dot = hat(b2_c_dot) * b3_c + hat(b2_c)*b3_c_dot;
                b1_c_ddot = hat(b2_c_ddot) * b3_c + 2 * hat(b2_c_dot) * b3_c_dot + hat(b2_c) * b3_c_ddot;
                R_c = [b1_c,b2_c,b3_c];
                R_c_dot = [b1_c_dot, b2_c_dot, b3_c_dot];
                R_c_ddot = [b1_c_ddot, b2_c_ddot, b3_c_ddot];
   
                W_c = vee(R_c' * R_c_dot);
                W_c_dot = vee(R_c' * R_c_ddot - hat(W_c)^2);

%                 W3 = dot(R_now * uav.e3, R_c * W_c);
%                 W3_dot = dot(R_now * uav.e3, R_c * W_c_dot) + dot(R_now * hat(W_now) * uav.e3, R_c * W_c);
                
                W_hat = hat(W_now);
                eR = 0.5*vee((R_c'*R_now-R_now'*R_c));     % error R
                eW = W_now-R_now'*R_c*W_c;

                
                
                if type == "origin"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot);
                elseif type == "EMK"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot) - [uav.pc_2_mc(2);-uav.pc_2_mc(1);0]*f;
                elseif type == "adaptive"
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot) - obj.theta*f;
                    theta_hat_dot = obj.gamma*f*(eW+obj.c2*eR);  
                    obj.theta = obj.theta + theta_hat_dot;
                    obj.theta(3) = 0;
                    
                    obj.theta = obj.theta + theta_hat_dot ;
                    obj.theta(3) = 0;
                    disp("theta")
                    disp(obj.theta);
                elseif type == "ICL"

                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot) - obj.theta*f;

                    
                    Y_omega = [                        0 , -W_now(2,1)*W_now(3,1) ,  W_now(2,1)*W_now(3,1) ;
                                   W_now(1,1)*W_now(3,1) ,                      0 , -W_now(1,1)*W_now(3,1) ;
                                  -W_now(1,1)*W_now(2,1) ,  W_now(1,1)*W_now(2,1) ,                      0 ;];
                        
                    Y_omega_J = Y_omega*[uav.J(1,1);uav.J(2,2);uav.J(3,3)];
                    obj.y       = obj.y + f*uav.dt - obj.Y_array(1)*uav.dt;
                    obj.y_omega = obj.y_omega + Y_omega_J*uav.dt - obj.Y_omega_array(:,1)*uav.dt + [(W_now(1)-obj.W_array(1,1))*uav.J(1,1);(W_now(2)-obj.W_array(2,1))*uav.J(2,2);(W_now(3)-obj.W_array(3,1))*uav.J(3,3)];
                    obj.M_hat   = obj.M_hat + M*uav.dt -obj.M_array(:,1)*uav.dt;

%                     Y_omega_J = Y_omega*[uav.J(1,1);uav.J(2,2);uav.J(3,3)];
%                     obj.y       = obj.y + f - obj.Y_array(1);
%                     obj.y_omega = obj.y_omega + Y_omega_J - obj.Y_omega_array(:,1) + [(W_now(1)-obj.W_array(1,1))*uav.J(1,1);(W_now(2)-obj.W_array(2,1))*uav.J(2,2);(W_now(3)-obj.W_array(3,1))*uav.J(3,3)];
%                     obj.M_hat   = obj.M_hat + M -obj.M_array(:,1);
                    
                    for i= 1:obj.integral_times_discrete-1
                        obj.Y_array(i) = obj.Y_array(i+1);
                        obj.Y_omega_array(:,i) = obj.Y_omega_array(:,i+1);
                        obj.M_array(:,i) = obj.M_array(:,i+1);
                        obj.W_array(:,i) = obj.W_array(:,i+1);
                        obj.R_array(:,i) = obj.R_array(:,i+1);
                    end
                    obj.Y_array(obj.integral_times_discrete)         = f;
                    obj.Y_omega_array(:,obj.integral_times_discrete) = Y_omega_J;
                    obj.M_array(:,obj.integral_times_discrete)       = M;
                    obj.W_array(:,obj.integral_times_discrete)       = W_now;
                    obj.R_array(:,obj.integral_times_discrete)       = R_now;
%                     disp("obj.m array");
%                     disp(obj.M_array);
%                         disp(" obj.y_omega");
%                         disp( obj.y_omega);
%                         disp("obj.y ");
%                         disp(obj.y);
%                         disp(" obj.M_hat");
%                         disp(obj.M_hat);
                    if iteration > obj.integral_times_discrete
                        for i= 1:obj.N-1
                            obj.sigma_M_hat_array(:,i) = obj.sigma_M_hat_array(:,i+1);
                            obj.sigma_y_omega_array(:,i) = obj.sigma_y_omega_array(:,i+1);
                            obj.sigma_y_array(i) = obj.sigma_y_array(i+1);
                        end
                        
                        obj.sigma_M_hat_array(:,obj.N) = obj.M_hat;
                        obj.sigma_y_omega_array(:,obj.N) = obj.y_omega;
                        obj.sigma_y_array(obj.N) = obj.y;
                        
                        x = zeros(3,1);

                        for i=1:obj.N
                                x = x + obj.sigma_y_array(i)*( obj.sigma_y_omega_array(:,i) - obj.sigma_M_hat_array(:,i) - obj.sigma_y_array(i)*obj.theta );
%                                 sigma_y_omega = sigma_y_omega +  obj.sigma_y_omega_array(:,i);
%                                 sigma_M = sigma_M - obj.sigma_M_hat_array(:,i);
%                                 sigma_y_theta = sigma_y_theta + obj.sigma_y_array(i)*obj.theta;
                                
                        end

%                         disp("x");
%                         disp(x);
                        theta_hat_dot = obj.gamma*f*(eW+obj.c2*eR) + obj.gamma*obj.k_icl * x;

                    
                    else
                        theta_hat_dot = obj.gamma*f*(eW+obj.c2*eR);
                    end
                    obj.theta = obj.theta + theta_hat_dot ;
                    obj.theta(3) = 0;
                    %disp("theta")
                    %disp(obj.theta);
                end
%                 disp("M")
%                 disp(M);
                control(1) = f;
                control(2:4) = M;

              end
   end
end
