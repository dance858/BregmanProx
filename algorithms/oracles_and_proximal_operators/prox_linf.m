% Copyright Laurent Condat
% v1.0, May 11, 2018
function x = prox_linf (y, lambda)
	tmp = abs(y(:));
	if sum(tmp)<=lambda, x = zeros(size(y)); return; end
	tau = max((cumsum(sort(tmp,1,'descend'))-lambda)./(1:length(tmp))');
	x = max(min(y, tau), -tau);
end