function [G] = combiner_interference_ic(H, W)
	[N_e, K] = deal(size(W, 2), size(W, 4));
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
    [U, ~] = pageeig(sum(T, 4) - T(:, :, logical(eye(K))));
    G = pagectranspose(U(:, 1 : N_e, :));
end
