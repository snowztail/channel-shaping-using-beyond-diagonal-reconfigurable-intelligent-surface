function [Theta, H] = scatter_rate_ic(H_d, H_f, H_b, W, Theta, L, P_n, rho)
	G = length(Theta) / L;
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	J = rho' * rate_mimo_ic(H, W, P_n);
	while ~iter.converge
		[iter.J, iter.G_r, iter.D] = deal(J, G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) rho' * rate_mimo_ic(H_d + pagemtimes(pagemtimes(H_b(:, S, :, :), Theta_g), H_f(S, :, :, :)) + pagemtimes(pagemtimes(H_b(:, S_c, :, :), Theta(S_c, S_c)), H_f(S_c, :, :, :)), W, P_n);
			G_e(S, S) = gradient_euclidean(H, H_f(S, :, :, :), H_b(:, S, :, :), W, P_n, rho);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		J = rho' * rate_mimo_ic(H, W, P_n);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end

function [G_e] = gradient_euclidean(H, H_f, H_b, W, P_n, rho)
	[N_r, N_e, K] = deal(size(H, 1), size(W, 2), size(W, 4));
	T = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H, W), 'ctranspose');
	Q = sum(T, 4) - T(:, :, logical(eye(K))) + P_n * eye(N_r);

	F = pagemtimes(H(:, :, logical(eye(K))), W(:, :, 1 : end));
	E = pageinv(eye(N_e) + pagemtimes(pagemrdivide(pagectranspose(F), Q), F));

	A = pagemtimes(pagemtimes(H, W), 'none', pagemtimes(H_f, W), 'ctranspose');
	G_e = 0;
	for k = 1 : K
		G_e = G_e + rho(k) * H_b(:, :, k)' * inv(Q(:, :, k)) * H(:, :, k, k) * W(:, :, :, k) * E(:, :, k) * W(:, :, :, k)' * (H_f(:, :, :, k)' + H(:, :, k, k)' * inv(Q(:, :, k)) * sum(A(:, :, k, :), 4));
	end
end
