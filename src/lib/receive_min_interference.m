function [G] = receive_min_interference(H, W)
	[N_e, N_r, K] = deal(size(W, 2), size(H, 1), size(W, 3));
	G = zeros(N_e, N_r, K);
	Q = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	for k = 1 : K
		[U, ~] = eigs(sum(Q(:, :, k, 1 : end ~= k), 4), N_e, 'smallestabs');
		G(:, :, k) = U';
	end
end
