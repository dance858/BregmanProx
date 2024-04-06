function [value] = evaluate_regularizer(X)

value = 0;
m = size(X, 1);
p = size(X, 2)/m - 1;
a = zeros(1, 2*p+1);
% Gather entries
for i = 1:m-1
    for j = i+1:m
        a(1) = X(i, j);
        for k = 2:p+1
            a(k) = X(i, (k-1)*m + j);
            a(p+k) = X(j, (k-1)*m + i);
        end
        value = value + max(abs(a));
    end
end







end