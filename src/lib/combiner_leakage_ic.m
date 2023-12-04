function [G] = combiner_leakage_ic(H, W)
	[d, K] = deal(size(W, 2), size(W, 4));
	Q = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
    [U, ~] = pageeig(sum(Q, 4) - Q(:, :, logical(eye(K))));
    G = pagectranspose(U(:, 1 : d, :));
end
