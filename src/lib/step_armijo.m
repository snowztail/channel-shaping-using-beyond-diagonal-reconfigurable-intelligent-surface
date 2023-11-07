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
