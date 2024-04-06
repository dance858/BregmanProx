function [X_k, objs, stats_projs, opt_lambdas, pri_res, dual_res] = ... 
    BSP_sparse_spectral_est(X, X_grad_phi, Y, A, rhs, max_iter, tol_subp, ...
                            tau, sigma, gamma, verbose, ptol, dtol)
m = size(X, 1); X_k = X; X_k_minus_1 = X; X_k_grad_phi = X_grad_phi;
opt_lambda = 0;

% Quantities for tracking the progress
stats_projs = zeros(2, 0); opt_lambdas = []; objs = []; pri_res = []; dual_res = [];

for k = 1:max_iter
    % evaluate obj
    obj =  0.5*((norm(X_k(:, 1:m) - A(:, 1:m),'fro')^2)...
          + 2*norm(X_k(:, m+1:end) - A(:, m+1:end), 'fro')^2) ...
          + gamma*evaluate_regularizer(X_k);
    objs(end + 1) =  obj;
    
    if ~rem(k, 20) && verbose
        fprintf('iter/pri_res/dual_res %i / %e / %e\n', k, pri_res(end), dual_res(end));
    end
    
    % Step 1
    S_k = 2*X_k - X_k_minus_1;
    
    % Step 2
    Z = Y + sigma*S_k;
    Y_cand = Z - sigma*prox_sparse_euclidean(1/sigma*Z, 1/sigma, A, gamma); 
        
    % Step 3
    temp = tau*Y_cand - X_k_grad_phi;
    [X_cand, grad_phi_X_cand, ~, stats, opt_lambda] = ...
        MISproj(temp, rhs, tol_subp, opt_lambda);

    % Store some statistics
    stats_projs(:, end+1) = stats; opt_lambdas(end+1) = opt_lambda;

    % Track primal and dual residual
    p_res_vec = 1/tau*(X_k_grad_phi - grad_phi_X_cand);
    d_res_vec = 1/sigma*(Y - Y_cand) + (2*X_k - X_k_minus_1 - X_cand);
   
    pri_res(end+1) = sqrt(norm(p_res_vec(:, 1:m), 'fro')^2 + ...
                          2*norm(p_res_vec(:, m+1:end), 'fro')^2);
    dual_res(end+1) = sqrt(norm(d_res_vec(:, 1:m), 'fro')^2 + ...
                           2*norm(d_res_vec(:, m+1:end), 'fro')^2);
    
    % Update quantities for next iteration.
    X_k_minus_1 = X_k; X_k = X_cand; X_k_grad_phi = grad_phi_X_cand; Y = Y_cand;
    
    % Check termination criteria.
    if (pri_res(end) <= ptol && dual_res(end) <= dtol)
        break
    end
end
end