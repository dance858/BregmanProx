function [sol,grad_phi,phi,stats, opt_lambda] = MISproj(C, rhs, tol, lambda0)
% Solves the problem
% minimize      <C, Y> + phi(Y)
% subject to    Tr(Y0) = rhs.
%  Note: Y = [Y0, Y1, ..., Yp], C = [C0, C1, ..., Cp]
%        <C, Y> = Tr(C0*Y0) + 2*[Tr(C1'*Y1) + ... + Tr(Cp'*Yp)] etc.
%        
%  INPUT:
%        C: matrix of size m x m(p+1)
%      tol: tolerance for Newton's method (the algorithm terminates when 
%           the derivative of the dual objective has magnitude less than
%           tol)
%  lambda0: initial guess. 

maxiter = 30; m = size(C, 1); p = size(C, 2)/m - 1;
E = zeros(m, (p+1)*m); E(:, 1:m) = eye(m);
stats = zeros(2,1); %[num of iter, num of times fstchol is called]

% Compute initial guess for dual variable lambda, and a lower bound on 
% the optimal value of lambda. Also compute upper triangular U such that
% (Toep(C^T) + lambda*I)^{-1} = U*U'.
[lambda, lb_opt, ~, U, num_fstchol_called] = initial_guess(C, m, p, lambda0);
stats(2) = stats(2) + num_fstchol_called;

% Evaluate dual objective function. To evaluate the log-determinant in 
% a numerically safe way we use Cholesky.
s_val = rhs*lambda + 2*sum(log(diag(chol(U(1:m, :)*U(1:m, :)'))));


for ii = 1:maxiter
    if ii == maxiter - 1
       disp('Max iterations in MISproj!') 
    end
    % Compute first and second derivative of s in the current iterate.
    [s_prime, s_prime_prime] = s_diff(U, m, rhs);
   
    % Check termination criteria (Newton decrement)
    if abs(s_prime^2/(2*s_prime_prime)) < tol 
        break;
    end
                              
    % Take a Newton step. 
    dir = -s_prime/s_prime_prime;
    [step_size, U, s_val, lb_opt, num_fstchol_called] = ...
        compute_step_size(lambda, dir, lb_opt, s_prime, s_val, C, E, rhs);
    lambda = lambda + step_size*dir;
    stats(2) = stats(2) + num_fstchol_called;
end

% First recover the spectral factor B of the optimal Y, and then Y itself. 
% Note that R has in fact already been computed inside 'compute_s_derivatives',
% so possible to refactor. 
R = U*U(1:m, :)'; 
B = zeros(m*(p+1), m);
B(1:m, :) = chol(R(1:m, :))';
B(m:end, :) = R(m:end, :)/B(1:m, :)';
sol = block_diag_sum_two_transpose(B*B', m, p);

%B_test = spectral_factorization_via_DARE(sol);
% Recover the value and gradient of phi at the optimal Y.
% When evaluating phi we take into account that B0 is upper triangular.
grad_phi = -C;
grad_phi(:, 1:m) = grad_phi(:, 1:m) - lambda*eye(m);
phi = -2*sum(log(diag(B(1:m, :))));  
opt_lambda = lambda;
% Update statistics for solver.
stats(1) = ii;

end