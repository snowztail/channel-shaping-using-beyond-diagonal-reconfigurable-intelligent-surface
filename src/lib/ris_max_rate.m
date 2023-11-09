function [Theta, H] = ris_max_rate(H_d, H_f, H_b, Theta, G, Q, P_n)
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	R = rate_mimo(H, Q, P_n);
	while ~iter.converge
		[iter.R, iter.G_r, iter.D] = deal(R, G_r, D);
		for g = 1 : G
			S = (g - 1) / G * length(Theta) + 1 : g / G * length(Theta);
			fun = @(Theta) rate_mimo(channel_aggregate(H_d, H_f(S, :), H_b(:, S), Theta), Q, P_n);
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), Q, P_n);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		R = rate_mimo(H, Q, P_n);
		iter.converge = (abs(R - iter.R) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end

function [G_e] = gradient_euclidean(H, H_f, H_b, Q, P_n)
	G_e = H_b' / (eye(size(H, 1)) + H * (Q / P_n) * H') * H * (Q / P_n) * H_f';
end
