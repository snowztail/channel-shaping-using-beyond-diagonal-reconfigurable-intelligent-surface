function [G] = receive_min_interference(H, W)
	[N_e, N_r, K] = deal(size(W, 2), size(H, 1), size(W, 3));
	G = zeros(N_e, N_r, K);
	Q = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	for k = 1 : K
		[U, ~] = eigs(sum(Q(:, :, k, 1 : end ~= k), 4), N_e, 'smallestabs');
		G(:, :, k) = U';
	end

    for i = 1 : K
        M = 0;
        for j = setdiff(1 : K, i)
            M = M + H(:, :, i, j) * W(:, :, j) * W(:, :, j)' * H(:, :, i, j)';
        end
        [U, ~] = eigs(M, N_e, 'smallestabs');
        g(:, :, i) = U';
    end
end
