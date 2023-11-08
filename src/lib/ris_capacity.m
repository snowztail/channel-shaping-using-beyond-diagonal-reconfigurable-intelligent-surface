function [Theta, H] = ris_capacity(H_d, H_f, H_b, Theta, G, Q, P_n)
	fun = @(Theta) rate_mimo(channel_aggregate(H_d, H_f, H_b, Theta), Q, P_n);
	
	L = length(Theta) / G;
	[iter.converge, iter.counter, iter.R, iter.tolerance] = deal(false, 0, 0, 1e-6);
	while ~iter.converge
		for g = 1 : G
			S = (g-1)*L+1 : g*L;
			G_e = gradient_euclidean(H_d, H_f, H_b, Theta, S, Q, P_n);
			G_r = gradient_riemannian(Theta, G_e, S);
			D = direction_conjugate(G_r, S, iter);
			mu = step_armijo(fun, Theta, D, S);
			Theta = update_geodesic(Theta, D, mu, S);
			H = channel_aggregate(H_d, H_f, H_b, Theta);
			R = rate_mimo(H, Q, P_n);
		end
		iter.converge = (abs(R - iter.R) <= iter.tolerance);
		iter.counter = iter.counter + 1;
		[iter.G_r, iter.D, iter.Theta, iter.R] = deal(G_r, D, Theta, R);
	end
end

function [G_e] = gradient_euclidean(H_d, H_f, H_b, Theta, S, Q, P_n)
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	G_e(S, S) = H_b(:, S)' / (eye(size(H, 1)) + H * (Q / P_n) * H') * H * (Q / P_n) * H_f(S, :)';
end

function [G_r] = gradient_riemannian(Theta, G_e, S)
	G_r(S, S) = G_e(S, S) * Theta(S, S)' - Theta(S, S) * G_e(S, S)';
end

function [D] = direction_conjugate(G_r, S, iter)
	L = length(S);
	
	% * Periodic conjugate direction restart by steepest ascent
	if mod(iter.counter, L ^ 2) == 0
		gamma = 0;
		iter.D(S, S) = zeros(L);
	else
		gamma = trace((G_r(S, S) - iter.G_r(S, S)) * G_r(S, S)') / trace(iter.G_r(S, S) * iter.G_r(S, S)');
		
		% * If the conjugate direction is not ascent, restart
		if real(trace((G_r(S, S) + gamma * iter.D(S, S))' * G_r(S, S))) < 0
			gamma = 0;
		end
	end
	D(S, S) = G_r(S, S) + gamma * iter.D(S, S);
end

function [mu] = step_armijo(fun, Theta, D, S)
	O = fun(Theta);
	mu = 1;
	T = eye(size(Theta));
	T(S, S) = expm(mu * D(S, S));
	
	% * Undershoot, double the step size
	while (fun(T ^ 2 * Theta) - O) >= (mu * 0.5 * trace(D(S, S) * D(S, S)'))
		mu = mu * 2;
		T(S, S) = T(S, S) ^ 2;
	end
	
	% * Overshoot, halve the step size
	while (fun(T * Theta) - O) < (0.5 * mu * 0.5 * trace(D * D')) && (mu >= eps)
		mu = mu * 0.5;
		T(S, S) = expm(mu * D(S, S));
	end
end

function [Theta] = update_geodesic(Theta, D, mu, S)
	Theta(S, S) = expm(mu * D(S, S)) * Theta(S, S);
end
