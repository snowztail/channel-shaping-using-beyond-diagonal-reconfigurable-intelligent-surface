function [R] = rate_mimo_ic(H, W, P_n)
	[N_r, K] = deal(size(H, 1), size(H, 3));
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
	X = eye(N_r) + pagemtimes(pageinv(Q), T(:, :, logical(eye(K))));
	R = log(squeeze(prod(pageeig(X))));
	R = real(R);
end
