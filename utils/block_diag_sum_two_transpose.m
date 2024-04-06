function [output] = block_diag_sum_two_transpose(A, m, p)
% Represents the D-operator (that sums block below the diagonal).
output = zeros(m, m*(p+1));

for block_col = 0:p
        output(:, 1:(p+1)*m - block_col*m) = output(:, 1:(p+1)*m - block_col*m) ...
            + A(block_col*m + 1:end, block_col*m + 1: (block_col + 1)*m)';
end
end