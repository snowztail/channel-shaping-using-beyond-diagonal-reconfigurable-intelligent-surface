function [W] = precoder_rate_pc(H, P_t, P_n)
	[~, S, V] = svd(H, 'vector');
	P_n = P_n ./ S' .^ 2;
	P_s = waterfill(P_t, P_n(~isinf(P_n)));
	d = nnz(P_s);
	W = V(:, 1 : d) * diag(sqrt(P_s(1 : d)));
	% Q = W * W';
end
