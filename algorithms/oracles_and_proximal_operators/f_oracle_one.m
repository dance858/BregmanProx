function [val, grad] = f_oracle_one(X, A, status)
% Represents the function f(X) = 0.5*(||X_0 - A_0||_F^2 + 2*sum_{k=1}^p ||X_k - A_k||_F^2). 
X_minus_A = X-A; m = size(X, 1); 
val = norm(X_minus_A(:, 1:m), 'fro')^2 + 2*norm(X_minus_A(:, m+1:end), 'fro')^2;
val = 0.5*val;
if status == 1
    grad = -1;
else
    grad = X_minus_A;
end
end
