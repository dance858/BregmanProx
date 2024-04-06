function [result] = prox_sparse_euclidean(X, tau, A, gamma)
% Computes prox_{tau*g}(X) where g(Y) = 0.5*||Y - A||^2 + gamma*r(Y) and
% r(Y) = (group sparsity regularizer).
% INPUT:
%       X: (X0 ... Xp), size m x m(p+1)
%       A: (A0 ... Ap), size m x m(p+1)
%
% Returns a matrix of size m x m(p+1) representing prox_{tau*g}(X).     
% 
% NOTE: The norm that appears in the definition of g is induced
% by the inner product in the space of matrix coefficients, ie.
% ||Y - A||^2 = ||Y0 - A0||_F^2 + 2*||Y1 - A1||_F^2 + ...                 
%                               + 2*||Yp - Ap||_F^2.
m = size(X, 1); p = size(X, 2)/m - 1; result = zeros(m, m*(p+1));

% Gather entries
uij = zeros(1, 2*p+1);
aij = zeros(1, 2*p+1); 
for i = 1:m-1
    for j = i+1:m
        uij(1) = X(i, j);
        aij(1) = A(i, j);
        for k = 2:p+1
            uij(k) = X(i, (k-1)*m + j);
            uij(p+k) = X(j, (k-1)*m + i);
            
            aij(k) = A(i, (k-1)*m + j);
            aij(p+k) = A(j, (k-1)*m + i);
        end
        prox_this_pair = 1/(1+tau)*prox_linf (uij+tau*aij, gamma*tau/2);
        result(i, j) = prox_this_pair(1);
        result(j, i) = prox_this_pair(1);
        for k = 2:p+1
            result(i, (k-1)*m + j) = prox_this_pair(k);
            result(j, (k-1)*m + i) = prox_this_pair(p+k);
        end
    end
end
mask = repmat(eye(m), 1, p+1);
result = result + 1/(1+tau)*mask.*(X + tau*A);
end