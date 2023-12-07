function [W] = precoder_initial_ic(H, d, P_t)
	K = size(H, 3);
	[~, ~, V] = pagesvd(H(:, :, logical(eye(K))));
	W = sqrt(P_t / d) * permute(V(:, 1 : d, :), [1, 2, 4, 3]);
end
