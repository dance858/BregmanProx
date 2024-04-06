function [step_size, U, s_val_new, lb_opt, num_fstchol_called] = ...
    compute_step_size(lambda, dir, lb_opt, s_prime, s_val, C, E, rhs)
% Line search. 
% INPUT:
%       lambda: current iterate
%          dir: search direction
%       lb_opt: lower bound on the optimal value of lambda
%      s_prime: first derivative of dual objective s
%        s_val: dual objective s evaluated in lambda
%            C:
%            N:
%          rhs:
% OUTPUT:
%            U: upper triangular matrix that satisfies
%               (Toep(C^T) + lambda*Toep(N^T))^{-1} = U*U'.

C_size = size(C);
m = C_size(1);
dir_der = s_prime * dir;
step_size = 1;
num_fstchol_called = 0;

% Line search parameters.
beta = 0.5;
alpha = 0.01;

% Make sure that the step size is such that 
% lambda + step_size * direction > lb_opt.       
while (lambda + step_size*dir) <= lb_opt
    step_size = step_size/2;
end

% Armijo backtracking. 
while true 
    new_lambda = lambda + step_size*dir;
    pd_flag = true;
    try
        [~,~,U] = fstchol(C' + new_lambda*E');
    catch
       pd_flag = false; 
    end
    num_fstchol_called = num_fstchol_called + 1;
    
    % If Cholesky fails we should increase lb_opt.
    if ~pd_flag
        lb_opt = new_lambda;
    else  
        s_val_new = rhs*new_lambda + 2*sum(log(diag(chol(U(1:m, :)*U(1:m, :)'))));
        if (s_prime < 0) || (s_val_new <= (s_val + alpha*step_size*dir_der) + eps)
            break
        end
    end
    step_size = beta*step_size;
end
end