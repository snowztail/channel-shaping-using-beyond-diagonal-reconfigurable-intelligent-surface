clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 2 .^ (4 : 4 : 8), 4);
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.power, number.antenna, number.realization] = deal(length(transmit.power), length(reflect.antenna), 1e2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for p = 1 : number.power
		transmit.beamformer.direct = precoder_rate(channel.direct, transmit.power(p), receive.noise);
		receive.rate.direct(p, r) = rate_mimo(channel.direct, transmit.beamformer.direct, receive.noise);
	end
	% * Have RIS
	for a = 1 : number.antenna
		reflect.bond = reflect.antenna(a);
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		clear scatter_rate;
		for p = 1 : number.power
			transmit.beamformer.asymmetric = precoder_rate(channel.direct, transmit.power(p), receive.noise);
			[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.asymmetric, receive.noise));
			while ~iter.converge
				[reflect.beamformer.asymmetric, channel.aggregate.asymmetric] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer.asymmetric, reflect.bond, receive.noise);
				transmit.beamformer.asymmetric = precoder_rate(channel.aggregate.asymmetric, transmit.power(p), receive.noise);
				receive.rate.aggregate.asymmetric(p, a, r) = rate_mimo(channel.aggregate.asymmetric, transmit.beamformer.asymmetric, receive.noise);
				iter.converge = (abs(receive.rate.aggregate.asymmetric(p, a, r) - iter.rate) / iter.rate <= iter.tolerance);
				iter.rate = receive.rate.aggregate.asymmetric(p, a, r);
				iter.counter = iter.counter + 1;
			end

			reflect.beamformer.symmetric.enforced = projection_symmetry(reflect.beamformer.asymmetric);
			channel.aggregate.symmetric.enforced = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer.symmetric.enforced);
			transmit.beamformer.symmetric.enforced = precoder_rate(channel.aggregate.symmetric.enforced, transmit.power(p), receive.noise);
			receive.rate.aggregate.symmetric.enforced(p, a, r) = rate_mimo(channel.aggregate.symmetric.enforced, transmit.beamformer.symmetric.enforced, receive.noise);

			transmit.beamformer.symmetric.projected = precoder_rate(channel.direct, transmit.power(p), receive.noise);
			[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.symmetric.projected, receive.noise));
			while ~iter.converge && iter.counter <= 5e2
				[reflect.beamformer.symmetric.projected, channel.aggregate.symmetric.projected] = scatter_rate_symmetric_projected(channel.direct, channel.forward, channel.backward, transmit.beamformer.symmetric.projected, reflect.bond, receive.noise);
				transmit.beamformer.symmetric.projected = precoder_rate(channel.aggregate.symmetric.projected, transmit.power(p), receive.noise);
				receive.rate.aggregate.symmetric.projected(p, a, r) = rate_mimo(channel.aggregate.symmetric.projected, transmit.beamformer.symmetric.projected, receive.noise);
				iter.converge = (abs(receive.rate.aggregate.symmetric.projected(p, a, r) - iter.rate) / iter.rate <= iter.tolerance);
				iter.rate = receive.rate.aggregate.symmetric.projected(p, a, r);
				iter.counter = iter.counter + 1;
			end
		end
	end
end
receive.rate.direct = mean(receive.rate.direct, ndims(receive.rate.direct));
receive.rate.aggregate.asymmetric = mean(receive.rate.aggregate.asymmetric, ndims(receive.rate.aggregate.asymmetric));
receive.rate.aggregate.symmetric.enforced = mean(receive.rate.aggregate.symmetric.enforced, ndims(receive.rate.aggregate.symmetric.enforced));
receive.rate.aggregate.symmetric.projected = mean(receive.rate.aggregate.symmetric.projected, ndims(receive.rate.aggregate.symmetric.projected));
save('data/pc_rate_symmetry.mat');

figure('Name', 'Achievable Rate vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
handle.rate.direct = plot(pow2db(transmit.power), receive.rate.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for a = 1 : number.antenna
	handle.rate.aggregate(1, a) = plot(pow2db(transmit.power), receive.rate.aggregate.asymmetric(:, a) / log(2), 'DisplayName', 'Asymmetric BD: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
	handle.rate.aggregate(2, a) = plot(pow2db(transmit.power), receive.rate.aggregate.symmetric.enforced(:, a) / log(2), 'DisplayName', 'Symmetric-Enforced BD: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
	handle.rate.aggregate(3, a) = plot(pow2db(transmit.power), receive.rate.aggregate.symmetric.projected(:, a) / log(2), 'DisplayName', 'Symmetric-Projected BD: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
end
style_plot(handle.rate.aggregate, 3);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_symmetry.fig');
matlab2tikz('../assets/simulation/pc_rate_symmetry.tex', 'width', '12cm', 'height', '9cm');


function [Theta, H] = scatter_rate_symmetric_projected(H_d, H_f, H_b, W, L, P_n)
	Theta = eye(size(H_f, 1));

	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[G_e, G_r, D] = deal(zeros(size(Theta)));
	[iter.converge, iter.tolerance, iter.counter, iter.R] = deal(false, 1e-3, 0, rate_mimo(H, W, P_n));
	while ~iter.converge
		[iter.G_r, iter.D] = deal(G_r, D);
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			S_c = setdiff(1 : length(Theta), S);
			fun = @(Theta_g) rate_mimo(H_d + H_b(:, S) * Theta_g * H_f(S, :) + H_b(:, S_c) * Theta(S_c, S_c) * H_f(S_c, :), W, P_n);
			G_e(S, S) = gradient_euclidean(H, H_f(S, :), H_b(:, S), W, P_n);
			G_r(S, S) = gradient_riemannian(Theta(S, S), G_e(S, S));
			D(S, S) = direction_conjugate(G_r(S, S), struct('G_r', iter.G_r(S, S), 'D', iter.D(S, S), 'counter', iter.counter));
			Theta(S, S) = projection_symmetry(step_armijo(fun, Theta(S, S), D(S, S)));
			H = channel_aggregate(H_d, H_f, H_b, Theta);
		end
		R = rate_mimo(H, W, P_n);
		iter.converge = (abs(R - iter.R) / iter.R <= iter.tolerance);
		iter.R = R;
		iter.counter = iter.counter + 1;
	end
end

function [G_e] = gradient_euclidean(H, H_f, H_b, W, P_n)
	G_e = H_b' * H * W / (eye(size(W, 2)) + (H * W)' * (H * W) / P_n) * W' * H_f' / P_n;
end

function [Q] = projection_symmetry(X)
	Q = (X + X.') / 2;
end
