function [Theta, H] = update_ris(H_d, H_f, H_b, Theta, G, Q, P_n)
	[iter.converge, iter.tolerance, iter.counter, iter.R] = deal(false, 1e-6, 0, 0);
	while ~iter.converge
		G_r = gradient_riemannian(H_d, H_f, H_b, Theta, Q, P_n);
		D = direction_conjugate(G_r, iter);
		D = diagonal_block(D, G);
		mu = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D);
		Theta = expm(mu * D) * Theta;
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		R = rate_mimo(H, Q, P_n);

		iter.converge = (abs(R - iter.R) <= iter.tolerance);
		iter.counter = iter.counter + 1;
		[iter.G_r, iter.D, iter.Theta, iter.R] = deal(G_r, D, Theta, R);
	end
end

function [G_r] = gradient_riemannian(H_d, H_f, H_b, Theta, Q, P_n)
	H = channel_aggregate(H_d, H_f, H_b, Theta);

	% * Gradient on the Euclidean space
	G_e = (H_b') / (eye(size(H, 1)) + H * (Q / P_n) * H') * (H * (Q / P_n) * H_f');

	% * Gradient on the Riemannian space
	G_r = G_e * Theta' - Theta * G_e';
end

function [D] = direction_conjugate(G_r, iter)
	% * Periodic conjugate direction restart by steepest ascent
	if mod(iter.counter, numel(G_r)) == 0
		gamma = 0;
		iter.D = zeros(size(G_r));
	else
		gamma = trace((G_r - iter.G_r) * G_r') / trace(iter.G_r * iter.G_r');

		% * If the conjugate direction is not ascent, restart
		if real(trace((G_r + gamma * iter.D)' * G_r)) < 0
			gamma = 0;
		end
	end
	D = G_r + gamma * iter.D;
end

function [B] = diagonal_block(A, n)
	B = zeros(size(A));
	for i = 1 : n : length(B)
		B(i : i + n - 1, i : i + n - 1) = A(i : i + n - 1, i : i + n - 1);
	end
end

function [mu] = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D)
	R = rate_mimo(channel_aggregate(H_d, H_f, H_b, Theta), Q, P_n);

	iter.mu = 1;
	iter.T = expm(iter.mu * D);

	% * Undershoot, double the step size
	while (rate_mimo(channel_aggregate(H_d, H_f, H_b, iter.T ^ 2 * Theta), Q, P_n) - R) >= (iter.mu * 0.5 * trace(D * D'))
		iter.mu = iter.mu * 2;
		iter.T = iter.T ^ 2;
	end

	% * Overshoot, halve the step size
	while (rate_mimo(channel_aggregate(H_d, H_f, H_b, iter.T * Theta), Q, P_n) - R) < (0.5 * iter.mu * 0.5 * trace(D * D')) && (iter.mu >= eps)
		iter.mu = iter.mu * 0.5;
		iter.T = expm(iter.mu * D);
	end

	mu = iter.mu;
end
