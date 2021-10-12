classdef z_axis_KF
    properties
        P;
        Q;
        R;
        r_mp_z;
        states = zeros(7,1);
        A = eye(7)
        y_telta;
    end
    methods
        function obj = KF(obj,r_p,v_p,r_pm_x,r_pm_y,F,m,x,v,R,omega,dt)
         r_mp_x = -r_pm_x;
         r_mp_y = -r_pm_y;
         R_now = reshape(R,3,3);
         f = F(1);

         % predict
         B = [obj.states(4:6,1);
             9.81*[0;0;1] - (f/m)*R_now*[0;0;1];
             0]*dt;
         F_matrix = eye(7);
         obj.states = obj.states + B;
         obj.P = F_matrix*obj.P*F_matrix' + obj.Q;
         
         % correct
         D = [  r_mp_x*R_now(1,1) + r_mp_y*R_now(1,2);
                r_mp_x*R_now(2,1) + r_mp_y*R_now(2,2);
                r_mp_x*R_now(3,1) + r_mp_y*R_now(3,2);
                - R_now(1,1)*omega(3)*r_mp_y + R_now(1,2)*omega(3)*r_mp_x + R_now(1,3)*(-omega(2)*r_mp_x + omega(1)*r_mp_y);
                - R_now(2,1)*omega(3)*r_mp_y + R_now(2,2)*omega(3)*r_mp_x + R_now(2,3)*(-omega(2)*r_mp_x + omega(1)*r_mp_y);
                - R_now(3,1)*omega(3)*r_mp_y + R_now(3,2)*omega(3)*r_mp_x + R_now(3,3)*(-omega(2)*r_mp_x + omega(1)*r_mp_y);
                ];
         H = [  1,0,0,0,0,0,R_now(1,3);
                0,1,0,0,0,0,R_now(2,3);
                0,0,1,0,0,0,R_now(3,3);
                0,0,0,1,0,0,R_now(1,1)*omega(2)-R_now(1,2)*omega(1);
                0,0,0,0,1,0,R_now(2,1)*omega(2)-R_now(2,2)*omega(1);
                0,0,0,0,0,1,R_now(3,1)*omega(2)-R_now(3,2)*omega(1);];
       obj.y_telta =  [r_p';v_p'] - D - H*obj.states;
            
            S = H*obj.P*H' + obj.R;
            K = obj.P*H'/S;
            obj.states = obj.states + K*obj.y_telta;
            obj.P = (eye(7)-K*H)*obj.P;
        disp('y telta');
        disp(obj.y_telta);
         disp('states');
         disp(obj.states)
        end
    end
end