function [R] = rate_ic(H, W, G, P_n)
	[N_r, L, K] = deal(size(H, 1), size(W, 2), size(W, 4));
	for k = 1 : K
		for l = 1 : L
			B = eye(N_r) - H(:, :, k, k) * W(:, l, 1, k) * W(:, l, 1, k)' * H(:, :, k, k)' / (P_n * L);
			for j = 1 : K
				for d = 1 : L
					B = B + H(:, :, k, j) * W(:, d, 1, j) * W(:, d, 1, j)' * H(:, :, k, j)' / (P_n * L);
				end
			end
			sinr(k, l) = G(l, :, k) * H(:, :, k, k) * W(:, l, 1, k) * W(:, l, 1, k)' * H(:, :, k, k)' * G(l, :, k)' / (G(l, :, k) * B * G(l, :, k)')  / (P_n * L);
			R(k, l) = real(log(1 + sinr(k, l)));
		end
	end
end
