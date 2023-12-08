function [G] = combiner_initial_ic(H, N_e)
	K = size(H, 3);
	[U, ~, ~] = pagesvd(H(:, :, logical(eye(K))));
	G = pagectranspose(U(:, 1 : N_e, :));
end
