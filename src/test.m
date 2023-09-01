clear; setup; config_indoor;

iterator = create_iterator;
while ~iterator.converged
	[channel.aggregate] = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
	[bs.covariance] = update_bs(channel.aggregate, power.transmit, power.noise);
	[ris.scatter] = update_ris(channel.direct, channel.forward, channel.backward, bs.covariance, power.noise);
end

function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	H = H_d + H_b * Theta * H_f;
end

function [R] = rate_mimo(H, Q, P_n)
	R = log(det(eye(size(H, 1)) + (H * Q * H') / P_n));
	R = real(R);
end

function [Q] = update_bs(H, P_t, P_n)
	[~, S, V] = svd(H);
	lambda = diag(S .^ 2)'; lambda(lambda <= eps) = [];
	P_s = waterfill(P_t, P_n ./ lambda);
	Q = V(:, 1 : length(lambda)) * diag(P_s) * V(:, 1 : length(lambda))';
	Q = 0.5 * (Q + Q');
end

function [Theta] = update_ris(H_d, H_f, H_b, Q, P_n)
	persistent init;
	if isempty(init)
		init.Theta = eye(size(H_f, 1));
	end
	Theta = init.Theta;
	
	iter = create_iterator;
	while ~iter.converged && ~iter.timeout
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		G_e = (H_b') / (eye(size(H, 1)) + H * (Q / P_n) * H') * (H * (Q / P_n) * H_f');
		G_r = G_e * Theta' - Theta * G_e';
		
		if mod(iter.counter, numel(Theta)) == 0
			D = G_r;
		else
			gamma = trace((G_r - prev.G_r) * G_r') / trace(prev.G_r * prev.G_r');
			D = G_r + gamma * prev.D;
			if real(trace(D' * G)) < 0
				D = G_r;
			end
		end
		
		iter.converged = (0.5 * trace(G_r * G_r')) <= iter.minGap;
		if iter.converged
			break;
		end
		
		mu = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D);
		Theta = expm(-mu * D) * Theta;
		
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		R = rate_mimo(H, Q, P_n);
	end
	
	init.Theta = Theta;
end

function [mu] = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D)
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	R = rate_mimo(H, Q, P_n);
	
	test.mu = -1;
	test.T = expm(-test.mu * D);
	
	test.H = channel_aggregate(H_d, H_f, H_b, test.T ^ 2 * Theta);
	test.R = rate_mimo(test.H, Q, P_n);
	while (R - test.R) <= (test.mu * 0.5 * trace(D * D'))
		test.mu = test.mu * 2;
		% test.T = expm(-test.mu * D);
		test.T = test.T ^ 2;
		
		test.H = channel_aggregate(H_d, H_f, H_b, test.T ^ 2 * Theta);
		test.R = rate_mimo(test.H, Q, P_n);
	end
	
	test.H = channel_aggregate(H_d, H_f, H_b, test.T * Theta);
	test.R = rate_mimo(test.H, Q, P_n);
	while (R - test.R) > (0.5 * test.mu * 0.5 * trace(D * D'))
		test.mu = test.mu * 0.5;
		test.T = expm(-test.mu * D);
		% test.T = test.T ^ 0.5;
		
		test.H = channel_aggregate(H_d, H_f, H_b, test.T * Theta);
		test.R = rate_mimo(test.H, Q, P_n);
	end
	
	mu = test.mu;
end

function [iterator] = create_iterator
	[iterator.converged, iterator.timeout] = deal(false);
	[iterator.minGap, iterator.maxRun] = deal(1e-4, 1e2);
	[iterator.reference, iterator.counter] = deal(0);
end
