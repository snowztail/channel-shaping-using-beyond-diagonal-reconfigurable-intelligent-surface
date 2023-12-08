function [W] = precoder_rate_ic(H, W, P_t, P_n, rho)
	[N_r, N_t, N_e, K] = deal(size(H, 1), size(H, 2), size(W, 2), size(W, 4));
	[iter.converge, iter.tolerance, iter.counter] = deal(false, eps, 0);
	iter.J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
		Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
		F = pagemtimes(H(:, :, logical(eye(K))), W(:, :, 1 : end));
		G = pagemrdivide(pagectranspose(F), Q + T(:, :, logical(eye(K))));
		V = shiftdim(rho, -2) .* (eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));
		Y = pagemtimes(pagemtimes(pagemtimes(G, H), 'ctranspose', V, 'none'), pagemtimes(G, H));
		Z = pagemtimes(pagemtimes(W, 'ctranspose', Y, 'none'), W);
		lambda = (shiftdim(pagetrace(P_n * pagemtimes(pagemtimes(G, 'ctranspose', V, 'none'), G)) + pagetrace(sum(Z, 4)), -1) - pagetrace(sum(Z, 3))) / P_t;
		X(:, :, 1, :) = pagemtimes(pagemtimes(G, H(:, :, logical(eye(K)))), 'ctranspose', V, 'none');
		W = pagemldivide(sum(Y, 3) + lambda .* eye(N_t), X);
		J = rho' * rate_mimo_ic(H, W, P_n);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
end

function [b] = pagetrace(A)
	b = sum(pageeig(A), 1);
end
