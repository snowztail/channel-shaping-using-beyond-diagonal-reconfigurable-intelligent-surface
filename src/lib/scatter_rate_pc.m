function [Theta, H] = scatter_rate_pc(H_d, H_f, H_b, W, Theta, L, P_n)
	G = length(Theta) / L;
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	R = rate_mimo_pc(H, W, P_n);
	while ~iter.converge
		[iter.R, iter.G_r, iter.D] = deal(R, G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) rate_mimo_pc(H_d + H_b(:, S) * Theta_g * H_f(S, :) + H_b(:, S_c) * Theta(S_c, S_c) * H_f(S_c, :), W, P_n);
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), W, P_n);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		R = rate_mimo_pc(H, W, P_n);
		iter.converge = (abs(R - iter.R) / iter.R <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end

function [G_e] = gradient_euclidean(H, H_f, H_b, W, P_n)
	G_e = H_b' * H * W / (eye(size(W, 2)) + (H * W)' * (H * W) / P_n) * W' * H_f' / P_n;
end
