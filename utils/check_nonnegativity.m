function [nonnegative] = check_nonnegativity(Y)
% Tests if Y belongs to the cone of non-negative matrix polynomials.
m = size(Y, 1);
p = size(Y, 2)/m - 1;
%cvx_solver mosek
cvx_solver SDPT3
cvx_begin sdp
    variable XCVX((p+1)*m, (p+1)*m) symmetric 
    XCVX >= 0;
  
    for block_row = 0:p
        temp = zeros(m, m);
        for block_col = 0:p-block_row
            temp = temp +  XCVX((block_row + block_col)*m + 1: (block_row + block_col + 1)*m, ...
                              block_col*m + 1: (block_col + 1)*m);
        end
        Y(:, block_row*m+1:(block_row+1)*m) == temp;
    end
cvx_end

if strcmp(cvx_status, 'solved')
    nonnegative = true;
else
    nonnegative = false;
end
end