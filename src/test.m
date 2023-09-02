clear; close; setup;

[nTxs, nSxs, nRxs] = deal(4, 16, 2);
[power.transmit, power.noise] = deal(db2pow(-20), db2pow(-100));
[ris.full, ris.single] = deal(create_ris(nSxs, 1), create_ris(nSxs, nSxs));
channel = create_channel_indoor(nTxs, nRxs, nSxs);
channel.full = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.full.scatter);
channel.single = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.single.scatter);

[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-3, 0, 0);
while ~iter.converge
	bs.full.covariance = update_bs(channel.full, power.transmit, power.noise);
	ris.full.scatter = update_ris(channel.direct, channel.forward, channel.backward, bs.full.covariance, power.noise);
	channel.full = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.full.scatter);
	rate.full = rate_mimo(channel.full, bs.full.covariance, power.noise);
	
	iter.converge = (abs(rate.full - iter.rate) <= iter.tolerance);
	iter.counter = iter.counter + 1;
	iter.rate = rate.full;
end

function [H] = channel_aggregate(H_d, H_f, H_b, Theta)
	H = H_d + H_b * Theta * H_f;
end

function [R] = rate_mimo(H, Q, P_n)
	R = log(det(eye(size(H, 1)) + (H * Q * H') / P_n));
	R = real(R);
end

function [Q] = update_bs(H, P_t, P_n)
	% * Water-filling power allocation
	[~, S, V] = svd(H);
	lambda = diag(S .^ 2)'; lambda(lambda <= eps) = [];
	P_s = waterfill(P_t, P_n ./ lambda);
	
	% * Eigenmode transmission
	Q = V(:, 1 : length(lambda)) * diag(P_s) * V(:, 1 : length(lambda))';
	Q = 0.5 * (Q + Q');
end

function [Theta] = update_ris(H_d, H_f, H_b, Q, P_n)
	persistent iter;
	if isempty(iter)
		iter.Theta = eye(size(H_f, 1));
	end
	Theta = iter.Theta;
	
	[iter.converge, iter.tolerance, iter.counter, iter.R] = deal(false, 1e-3, 0, 0);
	while ~iter.converge
		G_r = gradient_riemannian(H_d, H_f, H_b, Theta, Q, P_n);
		D = direction_conjugate(G_r, iter);
		mu = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D);
		Theta = expm(mu * D) * Theta;
		R = rate_mimo(channel_aggregate(H_d, H_f, H_b, Theta), Q, P_n);
		
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
	while (rate_mimo(channel_aggregate(H_d, H_f, H_b, iter.T * Theta), Q, P_n) - R) < (0.5 * iter.mu * 0.5 * trace(D * D'))
		iter.mu = iter.mu * 0.5;
		iter.T = expm(iter.mu * D);
	end
	
	mu = iter.mu;
end
