function [G] = combiner_initial_ic(H, d)
	K = size(H, 3);
	[U, ~, ~] = pagesvd(H(:, :, logical(eye(K))));
	G = pagectranspose(U(:, 1 : d, :));
end
