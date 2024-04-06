function [X, obj, stats_proj, all_t, objs, opt_lambdas, time, lower_bound] = ...
    closest_seq_via_FBPG_non_monotonic(A, rhs, m, p, max_iter, tol_subp, t, ...
                              tol)
% Uses a fast Bregman proximal gradient method to solve the problem 
%      minimize     ||Y0 - A0||_F^2 + 2*(||Y1-A1||_F^2 + ... +
%                   ||Yp - Ap||_F^2)
%      subject to   X \in K, Tr(X0) = rhs.
%
% INPUT:
%       A: size m x m(p+1)
% OUTPUT:
%      stats_all_projections: #Newton iter and # Cholesky solves per
%                              Bregman projection
%             all_step_sizes: accepted value of tau in 'line search'
%                   all_objs: the objective value per iteration
%            all_opt_lambdas: the optimal dual multiplier per Bregman projection
%                       time: solver time
                          
X = zeros(m, m*(p+1));                          % Initialize at feasible
X(:, 1:m) = rhs/m*eye(m);                       % X = (rhs/m*I, 0, ..., 0).
X_grad_phi = zeros(m, m*(p+1));                 % The gradient of phi at X is
X_grad_phi(:, 1:m) = -m/rhs*eye(m);             % (-m/rhs*I, 0, ..., 0).
X_val_phi = -m*log(rhs/m);                      % The value of phi at X is -m*log(rhs/m).
f_oracle = @(X, status) f_oracle_one(X, A, status);
[X, obj, stats_proj, all_t, objs, opt_lambdas, time, lower_bound] = ...
    FBPG_non_monotonic_step(X, X_grad_phi, X_val_phi, rhs, f_oracle, max_iter, ...
                            tol_subp, t, tol, A);
                   
end