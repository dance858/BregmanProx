clear; clc; 
addpath('../algorithms/')
addpath('../utils/')
addpath('../algorithms/oracles_and_proximal_operators/')
addpath('../algorithms/entropic_prox_operator/')

%% Settings for experiment and FBPG
num_of_instances = 3; m = 10; all_p = [20, 30, 40, 50, 60, 70, 80];
max_iter = 2000; tol_subp = 1e-8; accuracy = 1e-3; t = 1e-4;
L0 = 10; beta = 2; gamma = 1.2;

%% Solve instances with FBPG
for k = 1:length(all_p)
    p = all_p(k);
    for run = 1:num_of_instances
        % ---------------- 
        % Load data
        % ----------------
        load(sprintf("results/cvx_p=%d_m=%d_run=%d.mat", p, m, run), ...
             'spectral_estimate')
        rhs = trace(spectral_estimate(:, 1:m));
        
        % ----------------
        % Solve with FBPG variant
        % ----------------
        [X, obj, stats_proj, all_t, objs, opt_lambdas, time, lower_bound] = ...
        closest_seq_via_FBPG_non_monotonic(spectral_estimate, rhs, m, p, ...
        max_iter, tol_subp, t, accuracy);
    
        % ----------------
        % Solve with different FBPG variant. This variant essentially has
        % identical performance to the variant above.
        % ----------------
        [X_new, obj_new, stats_proj_new, all_L_new, objs_new,  ...
         opt_lambdas_new, time_new, lower_bound_new] = ...
        closest_seq_via_FBPG_new(spectral_estimate, rhs, m, p, max_iter, ...
        tol_subp, L0, accuracy, beta, gamma);
    
        % ----------------
        % Store data
        % ----------------
        FBPG_struct = struct('stats_proj', stats_proj, 'all_t', all_t, ...
               'objs', objs, 'lower_bound', lower_bound, 'opt_lambdas', ...
                opt_lambdas, 'time', time, 'm', m, 'p', p);
        FBPG_struct_new = struct('stats_proj', stats_proj_new, 'all_L', ...
            all_L_new, 'objs', objs_new, 'lower_bound', lower_bound_new, ...
            'opt_lambdas', opt_lambdas_new, 'time', time_new, 'm', m, 'p', p);
        save(sprintf("results/FBPG_p=%d_m=%d_run=%d.mat", p, m, run), ...
            'FBPG_struct', 'FBPG_struct_new')
    end
end