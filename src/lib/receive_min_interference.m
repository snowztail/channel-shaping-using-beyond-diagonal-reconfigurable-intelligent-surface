function [U] = receive_min_interference(H, V)
	[N_e, N_r, K] = deal(size(V, 2), size(H, 1), size(V, 3));
	U = zeros(N_e, N_r, K);
	Q = pagemtimes(pagemtimes(H, V), 'none', pagemtimes(H, V), 'ctranspose');
	for k = 1 : K
		[u, ~] = eigs(sum(Q(:, :, k, :), 4), N_e, 'smallestabs');
		U(:, :, k) = u';
	end
end
