function [W] = precoder_rate_ic(H, W, P_t, P_n, rho)
	% [N_r, N_t, N_e, K, P_t] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4), pagenorm(W, 'fro') .^ 2);
	[N_r, N_t, N_e, K] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4));
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		F = pagemtimes(H(:, :, logical(eye(K))), pageswap(W));
		G = pagemrdivide(pagectranspose(F), Q + T(:, :, logical(eye(K))));
		V = shiftdim(rho, -2) .* pagehermitize(eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));
		Y = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		Z = pagehermitize(pagemtimes(pagemtimes(W, 'ctranspose', Y, 'none'), W));
		lambda = (pageswap(pagetrace(P_n * pagehermitize(pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G))) + pagetrace(sum(Z, 4))) - pagetrace(sum(Z, 3))) ./ P_t;
		W = pagemldivide(sum(Y, 3) + lambda .* eye(N_t), pageswap(pagemtimes(pagemtimes(G, H(:, :, logical(eye(K)))), 'ctranspose', V, 'none')));
		J = rho' * rate_mimo_ic(H, W, P_n);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance) && all(abs(pagenorm(W, 'fro') .^ 2 - P_t) ./ P_t < iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
end

function [b] = pagetrace(A)
	b = sum(pageeig(A), 1);
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end

function [A] = pagehermitize(A)
	A = 0.5 * (A + pagectranspose(A));
end
