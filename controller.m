classdef controller
    properties
         kx=13;
         kv=10;
         kR = 10*eye(3);
         kW = 1*eye(3);
         %% adaptive
         theta = [0;0;0];
         gamma = 0.0005;
%         gamma = 0;
         c2 = 1
        %% ICL
        integral_times_discrete ;
        k_icl = 0.0001;
        N = 20;
        y_i;
        
        y;
        y_omega;
        M_hat;
        
        Y_array;
        Y_omega_array;
        M_array;
        W_array;

        sigma_y_array;
        sigma_y_omega_array;
        sigma_M_hat_array;
         
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
                    
                    obj.theta = obj.theta + theta_hat_dot ;
                    obj.theta(3) = 0;
                    disp("theta")
                    disp(obj.theta);
                elseif type == "ICL"
                    disp("obj.integral_times_discrete");
                    disp(obj.integral_times_discrete);
                    M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*desired_W - R_now'*R_c*desired_W_dot) - obj.theta*f;
                    disp("M");
                    disp(M);
                    
                    Y_omega = [                        0 , -W_now(2,1)*W_now(3,1) ,  W_now(2,1)*W_now(3,1) ;
                                   W_now(1,1)*W_now(3,1) ,                      0 , -W_now(1,1)*W_now(3,1) ;
                                  -W_now(1,1)*W_now(2,1) ,  W_now(1,1)*W_now(2,1) ,                      0 ;];
                        
                    Y_omega_J = Y_omega*[uav.J(1,1);uav.J(2,2);uav.J(3,3)];
                    obj.y       = obj.y + f - obj.Y_array(1);
                    obj.y_omega = obj.y_omega + Y_omega_J - obj.Y_omega_array(:,1) + [(W_now(1)-obj.W_array(1,1))*uav.J(1,1);(W_now(2)-obj.W_array(2,1))*uav.J(2,2);(W_now(3)-obj.W_array(3,1))*uav.J(3,3)];
                    obj.M_hat   = obj.M_hat + M -obj.M_array(:,1);

                    
                    for i= 1:obj.integral_times_discrete-1
                        obj.Y_array(i) = obj.Y_array(i+1);
                        obj.Y_omega_array(:,i) = obj.Y_omega_array(:,i+1);
                        obj.M_array(:,i) = obj.M_array(:,i+1);
                        obj.W_array(:,i) = obj.W_array(:,i+1);
                    end
                    obj.Y_array(obj.integral_times_discrete)         = f;
                    obj.Y_omega_array(:,obj.integral_times_discrete) = Y_omega_J;
                    obj.M_array(:,obj.integral_times_discrete)       = M;
                    obj.W_array(:,obj.integral_times_discrete)       = W_now;
                    disp("obj.m array");
                    disp(obj.M_array);
                        disp(" obj.y_omega");
                        disp( obj.y_omega);
                        disp("obj.y ");
                        disp(obj.y);
                        disp(" obj.M_hat");
                        disp(obj.M_hat);
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
                        sigma_y_omega = zeros(3,1);
                        sigma_M = zeros(3,1);
                        sigma_y_theta = zeros(3,1);
                        for i=1:obj.N
                                x = x + obj.sigma_y_array(i)*( obj.sigma_y_omega_array(:,i) - obj.sigma_M_hat_array(:,i) - obj.sigma_y_array(i)*obj.theta );
                                sigma_y_omega = sigma_y_omega +  obj.sigma_y_omega_array(:,i);
                                sigma_M = sigma_M - obj.sigma_M_hat_array(:,i);
                                sigma_y_theta = sigma_y_theta + obj.sigma_y_array(i)*obj.theta;
                                
                        end

                        disp("x");
                        disp(x);
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
