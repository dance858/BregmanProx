function [Tx] = block_toep(X)
% Given X = (X0, X1, ..., Xp) this function returns the block Toeplitz
% matrix
%
% T(X) = [X0  X1'   ...  Xp'
%         X1  X0      
%         
%         Xp  ...        X0]
% X should be represented as a cell array with blocks!

num_of_blocks = length(X);
non_symm_block_toep = cell2mat(X(toeplitz(1:num_of_blocks)));
Tx = tril(non_symm_block_toep) + tril(non_symm_block_toep)' - diag(diag(non_symm_block_toep));
end

% A = [1 0; 0 4]; B = [5, 6; 7 8]; C = [9 10; 11 12], Q = {A, B, C};
% [Tx] = block_toep(X)