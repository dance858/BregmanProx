function [output] = block_diag_sum_one(A, m, p)
% Represents the D-operator (that sums block below the diagonal).
output = zeros(m, m*(p+1));

%for block_row = 0:p
%    for block_col = 0:p-block_row
for block_col = 0:p
    for block_row = 0:p-block_col
        output(1:m, block_row*m+1:(block_row+1)*m) = ...
            output(1:m, block_row*m+1:(block_row+1)*m) + ...   
            A((block_row + block_col)*m + 1: (block_row + block_col + 1)*m, ...
                 block_col*m + 1: (block_col + 1)*m);
    end
        %output(:, block_row*p + 1: block_row*(p+1)) = 1; 
    
end

end