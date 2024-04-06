function [output] = prox_g_conjugate(tau, Z, C, gamma)
% Computes prox_{tau * g_conj)(Z) where g(X) = <C, X> + sparsity_reg.
output = Z - tau*prox_sparsity_reg(gamma/tau, 1/tau*(Z - C));
end