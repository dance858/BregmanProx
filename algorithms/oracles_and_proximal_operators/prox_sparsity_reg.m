function [output] = prox_sparsity_reg(gamma, X)
% Computes prox_r(A) where r(X) = gamma * (group sparsity regularizer)
m = size(X, 1);
p = size(X, 2)/m - 1;
YPROX = zeros(m, m*(p+1));
% Gather entries
for i = 1:m-1
    for j = i+1:m
        a = zeros(1, 2*p+1);
        a(1) = X(i, j);
        for k = 2:p+1
            a(k) = X(i, (k-1)*m + j);
            a(p+k) = X(j, (k-1)*m + i);
        end
        prox_this_pair = prox_linf(a, gamma/2);
        YPROX(i, j) = prox_this_pair(1);
        YPROX(j, i) = prox_this_pair(1);
        for k = 2:p+1
            YPROX(i, (k-1)*m + j) = prox_this_pair(k);
            YPROX(j, (k-1)*m + i) = prox_this_pair(p+k);
        end
    end
end
mask = repmat(eye(m), 1, p+1);
output = YPROX + mask.*X;
end