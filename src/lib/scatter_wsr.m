function [Theta, H] = scatter_wsr(H_d, H_f, H_b, W, L, P_n, rho)
	persistent iter;
	if isempty(iter)
		Theta = eye(size(H_f, 1));
	else
		Theta = iter.Theta;
	end

	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter, iter.J] = deal(false, 1e-3, 0, sum(rho .* rate_mimo(H, W, P_n), 3));
	while ~iter.converge && iter.counter <= 1e2
		[iter.G_r, iter.D] = deal(G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) sum(rho .* rate_mimo(H_d + pagemtimes(pagemtimes(H_b(:, S, :, :), Theta_g), H_f(S, :, :, :)) + pagemtimes(pagemtimes(H_b(:, S_c, :, :), Theta(S_c, S_c)), H_f(S_c, :, :, :)), W, P_n), 3);
			G_e(S, S) = gradient_euclidean(H, H_f(S, :, :, :), H_b(:, S, :, :), W, P_n, rho);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		J = sum(rho .* rate_mimo(H, W, P_n), 3);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
	iter.Theta = Theta;
end

function [G_e] = gradient_euclidean(H, H_f, H_b, W, P_n, rho)
	[N_r, N_e, K] = deal(size(H, 1), size(W, 2), size(W, 4));
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);
	F = pagemtimes(H(:, :, logical(eye(K))), pageswap(W));
	A = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H_f, W), 'ctranspose');
	G_e = sum(rho .* pagemtimes(pagemtimes(pagemrdivide(pagectranspose(H_b), Q), pagemrdivide(F, eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F))), pagemtimes(pageswap(pagectranspose(W)), pageswap(pagectranspose(H_f)) - pagemtimes(pagemrdivide(pagectranspose(H(:, :, logical(eye(K)))), Q), sum(A, 4) - A(:, :, logical(eye(K)))))), 3);
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end
