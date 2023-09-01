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
		R = rate_mimo(H, Q, P_n);
		
		G_e = (H_b') / (eye(size(H, 1)) + H * (Q / P_n) * H') * (H * (Q / P_n) * H_f');
		G_r = G_e * Theta' - Theta * G_e';
		D = G_r;
		
		% * Armijo step size
		mu = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D);
		
		
		% * Conjugate direction
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		R = rate_mimo(H, Q, P_n);
		
		G_e = (H_b') / (eye(size(H, 1)) + H * (Q / P_n) * H') * (H * (Q / P_n) * H_f');
		G_r = G_e * Theta' - Theta * G_e';
		% gamma = trace(G_r -)
	end
end

function [mu] = step_armijo(H_d, H_f, H_b, Theta, Q, P_n, D)
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	R = rate_mimo(H, Q, P_n);
	
	test.mu = -1;
	test.T = expm(-test.mu * D);
	test.T2 = test.T ^ 2;
	
	test.Theta = test.T2 * Theta;
	test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
	test.R = rate_mimo(test.H, Q, P_n);
	while (R - test.R) <= (test.mu * 0.5 * trace(D * D'))
		test.T = test.T2;
		test.T2 = test.T ^ 2;
		test.mu = test.mu * 2;
		test.Theta = test.T2 * Theta;
		test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
		test.R = rate_mimo(test.H, Q, P_n);
	end
	
	test.Theta = test.T * Theta;
	test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
	test.R = rate_mimo(test.H, Q, P_n);
	while (R - test.R) > (0.5 * test.mu * 0.5 * trace(D * D'))
		test.T = expm(-test.mu * D);
		test.mu = test.mu * 0.5;
		test.Theta = test.T * Theta;
		test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
		test.R = rate_mimo(test.H, Q, P_n);
	end
	
	mu = test.mu;
end


% function [mu] = step_armijo(H_d, H_f, H_b, Q, Theta, G_r, D, P_n)
% 	H = channel_aggregate(H_d, H_f, H_b, Theta);
% 	R = rate_mimo(H, Q, P_n);

% 	test.mu = -1;
% 	test.T = expm(-test.mu * D);
% 	test.T2 = test.T ^ 2;

% 	test.Theta = test.T2 * Theta;
% 	test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
% 	test.R = rate_mimo(test.H, Q);
% 	while (R - test.R) <= (test.mu * 0.5 * trace(D * D'))
% 		test.T = test.T2;
% 		test.T2 = test.T ^ 2;
% 		test.mu = test.mu * 2;
% 		test.Theta = test.T2 * Theta;
% 		test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
% 		test.R = rate_mimo(test.H, Q);
% 	end

% 	test.Theta = test.T * Theta;
% 	test.H = channel_aggregate(H_d, H_f, H_b, test.Theta);
% 	test.R = rate_mimo(test.H, Q);
% 	while (R - test.R) > (0.5 * test.mu * 0.5 * trace(D * D'))
% 		test.T = expm(-test.mu * D);
% 		test.mu = test.mu * 0.5;
% 	end

% 	flag = 1;
% end


function [iterator] = create_iterator
	[iterator.converged, iterator.timeout] = deal(false);
	[iterator.minGap, iterator.maxRun] = deal(1e-2, 1e2);
	[iterator.reference, iterator.counter] = deal(0);
end
