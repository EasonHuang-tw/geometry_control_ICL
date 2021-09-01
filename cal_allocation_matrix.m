function M = cal_allocation_matrix(d, c_tau)
    cos_45 = cosd(45);
    M = [     1,     1,      1,     1;
         -d*cos_45, d*cos_45, d*cos_45, -d*cos_45;
          d*cos_45, d*cos_45,-d*cos_45, -d*cos_45;
         -c_tau, c_tau, -c_tau, c_tau];
end