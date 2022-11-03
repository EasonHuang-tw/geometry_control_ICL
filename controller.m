classdef controller
    properties
         kx=13;
         kv=6;
         kR = 5*eye(3);
         kW = 1*eye(3);
         
          M = [0;0;0];
         %% adaptive
%           theta = [0.01;0.01;0.01;-0.0005;-0.002;0.001;0;0];
          theta = [0;0;0;0;0;0;0;0];
%           gamma =  diag([0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00001,0.00001])*0.1;

%          gamma =  diag([0.000008,0.000008,0.000005,0.000005,0.000005,0.000005,0.000005,0.000005]);
            gamma =  diag([0.000005,0.000005,0.000005,0.000005,0.000005,0.000005,0.000005,0.000005]);
%         gamma = 0;
         c2 = 6.5
        %% ICL
        last_M= [0;0;0];
        last_M2= [0;0;0];
        last_M3= [0;0;0];
        last_M4= [0;0;0];
        last_M5= [0;0;0];
        
        last_W= [0;0;0];
        last_W2= [0;0;0];
        last_W3= [0;0;0];
        last_W4= [0;0;0];
        last_W5= [0;0;0];
        last_f =0;
        last_f2 =0;
        last_f3 =0;
        last_f4 =0;
        last_f5 =0;
        
        last_R = [1 0 0;0 1 0;0 0 1]
%         k_icl =  diag([2160,2160,2160,2160,2160,2160,216,216])*100;
%         k_icl =  diag([2160,2160,2160,2160,2160,2160,216,216])*90;
         k_icl =  diag([2200,2220,2220,2220,2220,2220,2220,2220])*200/200;
%        k_icl = 506000;
        N = 20;      
        
        theta_hat_dot= [0;0;0;0;0;0;0;0];
        
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
                    obj.M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot);
                elseif type == "EMK"
                    obj.M = -obj.kR * eR - obj.kW*eW + cross(W_now,uav.J*W_now) - uav.J*(W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot) + f*[-uav.pc_2_mc(2);uav.pc_2_mc(1);0];
                elseif type == "adaptive"
                    o_b = W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot; %omega_bar
                    Y1 = [   0   -f  ;...
                            f   0   ;...
                            0   0   ];
                    theta2 = [uav.J(1,1);uav.J(2,2);uav.J(3,3);uav.J(1,2);uav.J(1,3);uav.J(2,3);uav.pc_2_mc(1);uav.pc_2_mc(2)];
                    Y2 = [ -o_b(1) ,-W_now(2)*W_now(3) ,W_now(2)*W_now(3) ,-o_b(2)-W_now(1)*W_now(3) ,-o_b(3)+W_now(1)*W_now(2), -W_now(3)^2 + W_now(2)^2 ;...
                           W_now(1)*W_now(3), -o_b(2) ,-W_now(1)*W_now(3) ,-o_b(1)+W_now(2)*W_now(3) ,-W_now(1)^2 + W_now(3)^2 , -o_b(3)-W_now(1)*W_now(2);...
                           -W_now(1)*W_now(2) ,W_now(1)*W_now(2), -o_b(3) ,-W_now(2)^2 + W_now(1)^2  ,-o_b(1)-W_now(2)*W_now(3), -o_b(2)+W_now(1)*W_now(3)];
                       Y = [Y2,Y1];
                       
                       obj.theta_hat_dot = -obj.gamma*Y'*(eW+obj.c2*eR);  
%                        disp(size( obj.theta_hat_dot));
                       obj.theta = obj.theta + obj.theta_hat_dot;

                       obj.M = -obj.kR * eR - obj.kW*eW + Y*obj.theta;
%                      result = Y*theta2- (cross(W_now,uav.J*W_now) - uav.J*o_b+ f*[-uav.pc_2_mc(2);uav.pc_2_mc(1);0])
%                         result = Y*obj.theta - Y*theta2
%                     disp("theta")
%                     disp(obj.theta);
                elseif type == "ICL"
                    theta2 = [uav.J(1,1);uav.J(2,2);uav.J(3,3);uav.J(1,2);uav.J(1,3);uav.J(2,3);uav.pc_2_mc(1);uav.pc_2_mc(2)];

                    o_b = W_hat*R_now'*R_c*W_c - R_now'*R_c*W_c_dot; %omega_bar
                    Y1 = [   0   -obj.last_f  ;...
                            obj.last_f   0   ;...
                            0   0   ];
                   Y2 = [ -o_b(1) ,-W_now(2)*W_now(3) ,W_now(2)*W_now(3) ,-o_b(2)-W_now(1)*W_now(3) ,-o_b(3)+W_now(1)*W_now(2), -W_now(3)^2 + W_now(2)^2 ;...
                           W_now(1)*W_now(3), -o_b(2) ,-W_now(1)*W_now(3) ,-o_b(1)+W_now(2)*W_now(3) ,-W_now(1)^2 + W_now(3)^2 , -o_b(3)-W_now(1)*W_now(2);...
                           -W_now(1)*W_now(2) ,W_now(1)*W_now(2), -o_b(3) ,-W_now(2)^2 + W_now(1)^2  ,-o_b(1)-W_now(2)*W_now(3), -o_b(2)+W_now(1)*W_now(3)];
                   Y = [Y2,Y1];
%                    Y*theta2 - (cross(W_now,uav.J*W_now) - uav.J*o_b+ obj.last_f*[-uav.pc_2_mc(2);uav.pc_2_mc(1);0])
                   Y_omega = [ 0                   ,-W_now(2)*W_now(3) ,W_now(2)*W_now(3)             ,-W_now(1)*W_now(3)      ,W_now(1)*W_now(2)          , -W_now(3)^2 + W_now(2)^2 ;...
                                W_now(1)*W_now(3)   , 0                 ,-W_now(1)*W_now(3)            , W_now(2)*W_now(3)      ,-W_now(1)^2 + W_now(3)^2   , -W_now(1)*W_now(2);...
                                -W_now(1)*W_now(2)  ,W_now(1)*W_now(2)  , 0                            ,-W_now(2)^2 + W_now(1)^2,-W_now(2)*W_now(3)         , W_now(1)*W_now(3)];
%                    W_dot = (W_now-obj.last_W);
                   W_dot = (W_now-obj.last_W5);
                   
                   W_dot_matrix = [W_dot(1)     ,0        ,0            ,W_dot(2) ,W_dot(3),0       ;...
                                       0        , W_dot(2),0            ,W_dot(1) ,0       ,W_dot(3);...
                                       0        , 0       ,    W_dot(3) ,0        ,W_dot(1),W_dot(2)];
                                   
%                     M_bar = obj.M*uav.dt;
                    M_bar = obj.last_M5*5*uav.dt;
                    y_W = Y_omega*5*uav.dt + W_dot_matrix;
                    Y1_icl = [   0   -obj.last_f5  ;...
                            obj.last_f5   0   ;...
                            0   0   ];
                    y_cl = [y_W,Y1_icl*5*uav.dt];
%                     y_cl = [y_W,Y1*uav.dt];
                    
%                     (cross(W_now,uav.J*W_now)*uav.dt+uav.J*W_dot) - y_cl(:,1:6)*theta2(1:6)
                    if iteration > obj.N
                        for i= 1:obj.N-1
                            obj.sigma_M_hat_array(:,i) = obj.sigma_M_hat_array(:,i+1);
                            obj.sigma_y_array(:,:,i) = obj.sigma_y_array(:,:,i+1);
                        end
                        
                        obj.sigma_M_hat_array(:,obj.N) = M_bar;
                        obj.sigma_y_array(:,:,obj.N) = y_cl;

                        x = zeros(8,1);
                        for i=2:obj.N
                                x = x + obj.sigma_y_array(:,:,i)'*(obj.sigma_M_hat_array(:,i) - obj.sigma_y_array(:,:,i)*obj.theta );
%                                 disp(obj.sigma_M_hat_array(:,i) - obj.sigma_y_array(:,:,i)*obj.theta)
                        end
                         obj.theta_hat_dot = -obj.gamma*Y'*(eW+obj.c2*eR) + obj.gamma*obj.k_icl * x;
                           %disp('theta hat dot')
                           %disp(theta_hat_dot)
                    else
                        for i= 1:obj.N-1
                            obj.sigma_M_hat_array(:,i) = obj.sigma_M_hat_array(:,i+1);
                            obj.sigma_y_array(:,:,i) = obj.sigma_y_array(:,:,i+1);
                        end
                        
                        obj.sigma_M_hat_array(:,obj.N) = M_bar;
                        obj.sigma_y_array(:,:,obj.N) = y_cl;
                        obj.theta_hat_dot = -obj.gamma*Y'*(eW+obj.c2*eR);  
                    end
                    obj.last_W5 = obj.last_W4;
                    obj.last_W4 = obj.last_W3;
                    obj.last_W3 = obj.last_W2;
                    obj.last_W2 = obj.last_W;
                    obj.last_W = W_now;
                    obj.last_f5 = obj.last_f4;
                    obj.last_f4 = obj.last_f3;
                    obj.last_f3 = obj.last_f2;
                    obj.last_f2 = obj.last_f;
                    obj.last_f = f;
                    obj.last_R = R_now;
                    obj.theta = obj.theta + obj.theta_hat_dot ;
                    obj.M = -obj.kR * eR - obj.kW*eW + Y*obj.theta;
                    
                    obj.last_M5 = obj.last_M4;
                    obj.last_M4 = obj.last_M3;
                    obj.last_M3 = obj.last_M2;
                    obj.last_M2 = obj.last_M;
                    obj.last_M = obj.M;
                    
                    
                end
%                 disp("M")
%                 disp(M);
                control(1) = f;
                control(2:4) = obj.M;

              end
   end
end
