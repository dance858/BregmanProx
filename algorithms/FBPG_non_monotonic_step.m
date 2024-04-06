function [X, f_val_X, stats_proj, all_t, objs, opt_lambdas, time, lower_bound] = ...
    FBPG_non_monotonic_step(X, X_grad_phi, X_val_phi, rhs, f_oracle, max_iter, ...
                            tol_subp, t, tol, A)
     
% A fast Bregman proximal gradient method for solving
%      minimize     f(X)
%      subject to   X \in K, Tr(X0) = rhs.
%
% INPUT:
%               X:  Starting point, size m x m(p+1)
%      X_grad_phi:  Gradient of entropy in starting point
%       X_val_phi:   Value of entropy in starting point
%        f_oracle:  Value and gradient oracle of f. The function call
%                   f_oracle(X, 1) should return the value of f at x, and
%                   f_oracle(X, 2) should return the value and gradient.                    
%        tol_subp: the tolerance used for evaluating the entropic
%                  prox operator. Typical value is 1e-7
%               t: initial estimate of 1/L. Typical value 1e-3 
%             tol: the algorithm terminates either when the 
%                   maximum number of iterations is reached or when 
%                   the duality gap is less than tol

tic;
% Initialization
m = size(X, 1); p = size(X, 2)/m - 1; theta = 1;
E = zeros(m, m*(p+1)); E(:, 1:m) = eye(m);
V = X; V_grad_phi = X_grad_phi; V_phi = X_val_phi;
[f_val_X, ~] = f_oracle(X, 1); best_UB = f_val_X;
BT = block_toep(mat2cell((X - A), m, m*ones(1, p+1))'); 
best_LB = -f_val_X + inner_prod_E(A, A - X, m) + rhs*min(eig(BT));

beta = 0.5; gamma = 1.2; % Line search parameters      
opt_lambda = 0;  % Initial guess for projection problem.
t_old = t/beta;

% Quantities for tracking the progress.
stats_proj = zeros(2, 0);
all_t = []; objs = [f_val_X];
opt_lambdas = []; 
lower_bound = [best_LB];

for k = 1:max_iter
    if  ~rem(k, 10)
        fprintf('iter/obj: %i / %f \n ', k,  objs(end))
    end
    
    t = gamma*t/beta;                  
    terminate_linesearch = false;
    while ~terminate_linesearch
       t = t*beta;
       theta = -theta^2*t/(2*t_old) + sqrt(theta^4*t^2/(4*t_old^2) + theta^2*t/t_old);
       Y = (1-theta)*X + theta*V;       
       [f_val_Y, f_grad_Y] = f_oracle(Y, 2);
       tau = t/theta;
       [V_cand, V_cand_grad_phi, V_cand_phi, stats, opt_lambda] = ...
          MISproj(tau*f_grad_Y - V_grad_phi, rhs, tol_subp, opt_lambda);   
       X_cand = (1-theta)*X + theta*V_cand;
       
       % Store some statistics.
       stats_proj(:, end+1) = stats; opt_lambdas(end+1) = opt_lambda;
       
       % Check descent criteria.
       [f_val_X_cand, ~] = f_oracle(X_cand, 1);
       bregman_dist_V_cand_and_V = V_cand_phi - V_phi - ...
                                   inner_prod_E(V_grad_phi, V_cand - V, m);
       
       rhs_descent_test = (1-theta)*f_val_X ...
            + theta*(f_val_Y + inner_prod_E(f_grad_Y, V_cand - Y, m)...
                     + 1/tau*bregman_dist_V_cand_and_V);
        
       if f_val_X_cand <= rhs_descent_test 
            terminate_linesearch = true;
            all_t(end+1) = t;
       end   
    end
    % Accept iterate
    X = X_cand;  V = V_cand; V_phi = V_cand_phi; V_grad_phi = V_cand_grad_phi;
    f_val_X = f_val_X_cand; best_UB = min(f_val_X, best_UB);  t_old = t;    
    objs(end+1) = f_val_X;
    
    % Evaluate new lower bound every 20th iteration
    if ~rem(k, 20)
        BT = block_toep(mat2cell((X - A), m, m*ones(1, p+1))'); 
        best_LB = max(-f_val_X + inner_prod_E(A, A - X, m) + rhs*min(eig(BT)), best_LB);
        lower_bound(end+1) = best_LB;
    else
        lower_bound(end + 1) = lower_bound(end); 
    end
    
    % Evaluate termination criteria
    if abs(best_UB - best_LB)/abs(best_UB) <= tol
       break
    end
end
time = toc;
end