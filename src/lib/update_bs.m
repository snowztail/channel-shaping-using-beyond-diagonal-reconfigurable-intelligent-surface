function [Q, P_s] = update_bs(H, P_t, P_n)
	% * Water-filling power allocation
	[~, S, V] = svd(H);
	lambda = diag(S .^ 2)'; lambda(lambda <= eps) = [];
	P_s = waterfill(P_t, P_n ./ lambda);
	
	% * Eigenmode transmission
	Q = V(:, 1 : length(lambda)) * diag(P_s) * V(:, 1 : length(lambda))';
	Q = 0.5 * (Q + Q');
end
