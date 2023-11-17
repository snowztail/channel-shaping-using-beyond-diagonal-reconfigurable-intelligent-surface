function [R] = rate_ic(H, W, G, P_n)
	[N_r, L, K] = deal(size(H, 1), size(W, 2), size(W, 4));
	R = zeros(K, L);
	for k = 1 : K
		for l = 1 : L
			T = pagemtimes(H(:, :, k, :), W);
			B = eye(N_r) + (sum(pagemtimes(T, 'none', T, 'ctranspose'), 4) - H(:, :, k, k) * W(:, l, 1, k) * W(:, l, 1, k)' * H(:, :, k, k)') / (P_n * L);
			sinr = G(l, :, k) * H(:, :, k, k) * W(:, l, 1, k) * W(:, l, 1, k)' * H(:, :, k, k)' * G(l, :, k)' / (G(l, :, k) * B * G(l, :, k)')  / (P_n * L);
			R(k, l) = real(log(1 + sinr));
		end
	end
end
