function [s_prime, s_prime_prime] = s_diff(U, m, rhs)
% Evaluates the first and second derivative of the dual objective s.
% Note: The function s is defined as
%       s(lambda) = log det(E'*[Toep(C^T) + lambda*I]^{-1} E) + rhs*lambda.
% 
% INPUT:
%       U: An upper triangular matrix U satisfying 
%         (Toep(C^T) + lambda*I)^{-1} = U*U';
%       m: dimension of the blocks

% Compute U*U'*E
temp1 = U*U(1:m, :)';

% Compute B = E'*U*U'*U*U'E.
B = temp1'*temp1;

% Factor E'*(Toep(C^T) + lambda*I)^{-1} E = R'*R
R = chol(U(1:m, :)*U(1:m, :)');

first_deriv_matrix = (R\(R'\B));

% Compute first derivative.
s_prime = rhs - trace(first_deriv_matrix);

% Compute U'*Toep(N^T)*U*U'*E
temp3 = U'*temp1;

% Compute E'*U*U'*Toep(N^T)*U*U'*Toep(N^T)*U*U'*E.
temp4 = temp3'*temp3;

s_prime_prime = 2*trace((R\(R'\temp4))) - trace(first_deriv_matrix*first_deriv_matrix);

end