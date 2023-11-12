function [Q] = transmit_max_rate(H, P_t, P_n)
	[~, S, V] = svd(H, 'vector');
    S(S <= 1e-6) = []; V = V(:, 1 : length(S));
	P_s = waterfill(P_t, P_n ./ S' .^ 2);
    Q = V * diag(P_s) * V';
end
