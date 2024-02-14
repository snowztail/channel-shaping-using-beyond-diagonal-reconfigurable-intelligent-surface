function [W] = precoder_rate(H, P_t, P_n)
	[~, S, V] = svd(H, 'vector');
	P_n = P_n ./ S' .^ 2;
	P_s = waterfill(P_t, P_n(~isinf(P_n)));
	N_e = nnz(P_s);
	W = V(:, 1 : N_e) * diag(sqrt(P_s(1 : N_e)));
end
