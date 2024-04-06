function [XCVX, YCVX, cvx_optval, cvx_cputime, solve_time, cvx_status] = ...
    closest_sequence_via_cvx(A, rhs, m, p, precision, solver)
% Calls CVX to solve the problem 
%      minimize    0.5*( ||Y0 - A0||_F^2 + 2*(||Y1-A1||_F^2 + ... +
%                   ||Yp - Ap||_F^2) )
%      subject to   X \in K, Tr(X0) = rhs.
%
% INPUT:
%       A: size m x m(p+1)
% OUTPUT:
%       XCVX: block Toeplitz matrix of size m(p+1) x m(p+1)
%       YCVX: the solution
%       cvx_optval, cvx_cputume, cvx_status: output from cvx
%       solve_time: measured yourself
if precision == "low"
    cvx_precision low
elseif precision == "high"
   cvx_precision high 
end

if solver == "sedumi"
    cvx_solver sedumi
elseif solver == "SDPT3"
    cvx_solver SDPT3
elseif solver == "scs"
    cvx_solver scs
end

cvx_begin sdp
    variable XCVX((p+1)*m, (p+1)*m) symmetric 
    variable YCVX(m, (p+1)*m)                  % Note how we define it
    XCVX >= 0;
    trace(YCVX(:, 1:m)) == rhs;
    obj =  square_pos(norm(YCVX(:, 1:m) - A(:, 1:m),'fro')) + 2 * square_pos(norm(YCVX(:, m+1:end) - A(:, m+1:end), 'fro'));
    obj = 0.5*obj;
    for block_row = 0:p
        temp = zeros(m, m);
        for block_col = 0:p-block_row
            temp = temp +  XCVX((block_row + block_col)*m + 1: (block_row + block_col + 1)*m, ...
                              block_col*m + 1: (block_col + 1)*m);
        end
        YCVX(:, block_row*m+1:(block_row+1)*m) == temp;
    end
    minimize(obj)
    tic;
cvx_end
solve_time = toc;
end