function [x] = projToSmplx(v)
% Input: v is n x 1 column matrix
% Output: x is the result of projection
% Description: project v into simplex
u = sort(v,'descend');
sv = cumsum(u);
ind = find(u > (sv - 1) ./ (1:length(u))', 1, 'last');
tau = (sv(ind) - 1) / ind;
x = max(v - tau, 0);
end