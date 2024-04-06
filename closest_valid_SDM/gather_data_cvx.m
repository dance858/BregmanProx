clear; clc;
addpath('../utils/')
store_data = false;

%%
num_of_instances = 3; m = 10; all_p = [20, 30, 40, 50, 60, 70, 80];
precision = "low";

for p = all_p
   for problem_instance = 1:num_of_instances
   % ---------------- 
   % Generate random spectral estimate
   % ----------------
   spectral_estimate = randn(m, m*(p+1));
   spectral_estimate(:, 1:m) = spectral_estimate(:, 1:m)'*spectral_estimate(:, 1:m);
   rhs = trace(spectral_estimate(:, 1:m));
   
   % ---------------- 
   % Solve with IPM
   % ----------------
   [XCVX, YCVX, cvx_optval, cvx_cputime, solve_time, cvx_status] = ...
       closest_sequence_via_cvx(spectral_estimate, rhs, m, p, precision, "SDPT3");
   
   % ---------------- 
   % Store data
   % ----------------
   cvx_struct = []; cvx_struct.cvx_optval = cvx_optval;
   cvx_struct.cvx_cputime = cvx_cputime; cvx_struct.solve_time = solve_time;
   cvx_struct.cvx_status = cvx_status; cvx_struct.m = m;
   cvx_struct.p = p;
   if store_data
        save(sprintf("results/cvx_p=%d_m=%d_run=%d.mat", p, m, problem_instance), ...
            'cvx_struct', 'spectral_estimate');
   end
   end
end


