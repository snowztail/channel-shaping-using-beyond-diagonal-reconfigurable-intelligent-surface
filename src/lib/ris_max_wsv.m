function [Theta, H] = ris_max_wsv(H_d, H_f, H_b, rho, Theta, G)
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-6, 0);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	J = rho' * svd(H);
	while ~iter.converge
		[iter.J, iter.G_r, iter.D] = deal(J, G_r, D);
		for g = 1 : G
			S = (g - 1) / G * length(Theta) + 1 : g / G * length(Theta);
			fun = @(Theta) rho' * svd(channel_aggregate(H_d, H_f(S, :), H_b(:, S), Theta));
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), rho);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		J = rho' * svd(H);
		iter.converge = (abs(J - iter.J) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
end

function [G_e] = gradient_euclidean(H, H_f, H_b, rho)
	[U, S, V] = svd(H);
	S(boolean(eye(length(rho)))) = rho;
	G_e = H_b' * U * S * V' * H_f';
end
