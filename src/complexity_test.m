clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, [64, 256], 4);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 4);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.weight = ones(channel.rank, 1);
[number.antenna, number.realization, flag.direct] = deal(length(reflect.antenna), 1e2, true);

for r = 1 : number.realization
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));

		tic;
		[objective.geodesic(a, r), counter.geodesic(a, r)] = scatter_singular_geodesic(channel.direct, channel.forward, channel.backward, channel.weight, reflect.bond);
		timer.geodesic(a, r) = toc;

		tic;
		[objective.nongeodesic(a, r), counter.nongeodesic(a, r)] = scatter_singular_nongeodesic(channel.direct, channel.forward, channel.backward, channel.weight, reflect.bond);
		timer.nongeodesic(a, r) = toc;
	end
end

[objective.geodesic, objective.nongeodesic] = deal(mean(objective.geodesic, ndims(objective.geodesic)), mean(objective.nongeodesic, ndims(objective.nongeodesic)));
[counter.geodesic, counter.nongeodesic] = deal(mean(counter.geodesic, ndims(counter.geodesic)), mean(counter.nongeodesic, ndims(counter.nongeodesic)));
[timer.geodesic, timer.nongeodesic] = deal(mean(timer.geodesic, ndims(timer.geodesic)), mean(timer.nongeodesic, ndims(timer.nongeodesic)));

save('data/complexity_test.mat');



function [G_e] = gradient_euclidean(H, H_f, H_b, rho)
	[U, S, V] = svd(H);
	S(logical(eye(length(rho)))) = rho;
	G_e = H_b' * U * S * V' * H_f';
end

function [J, r] = scatter_singular_geodesic(H_d, H_f, H_b, rho, L)
	Theta = eye(size(H_f, 1));
	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter, iter.J] = deal(false, 1e-4, 0, rho' * svd(H));
	while ~iter.converge
		[iter.G_r, iter.D] = deal(G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) rho' * svd(H_d + H_b(:, S) * Theta_g * H_f(S, :) + H_b(:, S_c) * Theta(S_c, S_c) * H_f(S_c, :));
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), rho);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo_geodesic(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		J = rho' * svd(H);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
	r = iter.counter;
end

function [J, r] = scatter_singular_nongeodesic(H_d, H_f, H_b, rho, L)
	Theta = eye(size(H_f, 1));
	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter, iter.J] = deal(false, 1e-4, 0, rho' * svd(H));
	while ~iter.converge
		[iter.G_r, iter.D] = deal(G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) rho' * svd(H_d + H_b(:, S) * Theta_g * H_f(S, :) + H_b(:, S_c) * Theta(S_c, S_c) * H_f(S_c, :));
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), rho);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = step_armijo_nongeodesic(fun, Theta(S, S), D(S, S));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		J = rho' * svd(H);
		iter.converge = (abs(J - iter.J) / iter.J <= iter.tolerance);
		iter.J = J;
		iter.counter = iter.counter + 1;
	end
	r = iter.counter;
end

function [Theta] = step_armijo_geodesic(fun, Theta, D)
	O = fun(Theta);
	mu = 1e2;
	T = expm(mu * D);
	% * Undershoot, double the step size
	while (fun(T ^ 2 * Theta) - O) >= (mu * 0.5 * trace(D * D'))
		mu = mu * 2;
		T = T ^ 2;
	end
	% * Overshoot, halve the step size
	while (fun(T * Theta) - O) < (0.5 * mu * 0.5 * trace(D * D')) && (mu >= eps)
		mu = mu * 0.5;
		T = expm(mu * D);
	end
	Theta = T * Theta;
end

function [Theta] = step_armijo_nongeodesic(fun, Theta, D)
	O = fun(Theta);
	mu = 1e2;
	W = (Theta + mu * D) * ((Theta + mu * D)' * (Theta + mu * D)) ^ -0.5;
	while (fun(W) - O) >= (mu * 0.5 * trace(D * D'))
		mu = mu * 2;
		W = (Theta + mu * D) * ((Theta + mu * D)' * (Theta + mu * D)) ^ -0.5;
	end
	% * Overshoot, halve the step size
	while (fun(W) - O) < (0.5 * mu * 0.5 * trace(D * D')) && (mu >= eps)
		mu = mu * 0.5;
		W = (Theta + mu * D) * ((Theta + mu * D)' * (Theta + mu * D)) ^ -0.5;
	end
	Theta = W;
end
