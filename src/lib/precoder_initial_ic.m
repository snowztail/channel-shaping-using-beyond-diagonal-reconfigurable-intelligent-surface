function [W] = precoder_initial_ic(H, N_e, P_t)
	K = size(H, 3);
	[~, ~, V] = pagesvd(H(:, :, logical(eye(K))));
	W = sqrt(P_t / N_e) * pageswap(V(:, 1 : N_e, :));
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end
