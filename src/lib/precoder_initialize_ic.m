function [W] = precoder_initialize_ic(H, d)
	K = size(H, 3);
	[~, ~, V] = pagesvd(H(:, :, logical(eye(K))));
	W = permute(V(:, 1 : d, :), [1, 2, 4, 3]);
end
