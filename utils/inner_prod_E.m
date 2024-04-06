function [out] = inner_prod_E(X, Y, m)
    out = trace(X(:, 1:m)*Y(:, 1:m)) + 2*trace(X(:, m+1:end)'*(Y(:, m+1:end)));
end