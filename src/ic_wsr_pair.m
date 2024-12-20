clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream] = deal(4, 64, 4, 2);
[network.coverage, network.pair] = deal(20, 2 : 2 : 10);
[transmit.power, receive.noise] = deal(db2pow(20), db2pow(-75));
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
[reflect.bond, channel.uncertainty] = deal([1, reflect.antenna], [1e-2, 1e-1, 5e-1]);
[number.bond, number.pair, number.uncertainty, number.realization] = deal(2, length(network.pair), length(channel.uncertainty), 1e2);

for r = 1 : number.realization
	for p = 1 : number.pair
		network.weight = ones(1, 1, network.pair(p));
		[network.coordinate.transmit, network.coordinate.reflect, network.coordinate.receive] = deal(distribution_disk(network.coverage, network.pair(p)), [0; 0], distribution_disk(network.coverage, network.pair(p)));
		channel.pathloss.direct = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.transmit), channel.pathloss.exponent.direct);
		channel.pathloss.forward = pageswap(pathloss(channel.pathloss.reference, vecnorm(network.coordinate.reflect - network.coordinate.transmit), channel.pathloss.exponent.forward));
		channel.pathloss.backward = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.reflect), channel.pathloss.exponent.backward);
		channel.direct.actual = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair(p), network.pair(p));
		channel.forward.actual = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna, transmit.antenna, 1, network.pair(p));
		channel.backward.actual = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna, network.pair(p), 1);
		% * No RIS
		transmit.beamformer = precoder_initialize(channel.direct.actual, transmit.stream, transmit.power);
		transmit.beamformer = precoder_wsr(channel.direct.actual, transmit.beamformer, transmit.power, receive.noise, network.weight);
		network.wsr.direct.estimate(p, r) = sum(network.weight .* rate_mimo(channel.direct.actual, transmit.beamformer, receive.noise), 3);
		network.wsr.direct.actual(p, r) = sum(network.weight .* rate_mimo(channel.direct.actual, transmit.beamformer, receive.noise), 3);
		% * Have RIS
		for u = 1 : number.uncertainty
			channel.forward.estimate = estimation_effect(channel.forward.actual, channel.uncertainty(u));
			channel.backward.estimate = estimation_effect(channel.backward.actual, channel.uncertainty(u));
			for b = 1 : number.bond
				clear scatter_wsr;
				transmit.beamformer = precoder_initialize(channel.direct.actual, transmit.stream, transmit.power);
				[iter.converge, iter.tolerance, iter.counter, iter.wsr] = deal(false, 1e-4, 0, sum(network.weight .* rate_mimo(channel.direct.actual, transmit.beamformer, receive.noise), 3));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate.estimate] = scatter_wsr(channel.direct.actual, channel.forward.estimate, channel.backward.estimate, transmit.beamformer, reflect.bond(b), receive.noise, network.weight);
					transmit.beamformer = precoder_wsr(channel.aggregate.estimate, transmit.beamformer, transmit.power, receive.noise, network.weight);
					network.wsr.aggregate.estimate(b, u, p, r) = sum(network.weight .* rate_mimo(channel.aggregate.estimate, transmit.beamformer, receive.noise), 3);
					iter.converge = (abs(network.wsr.aggregate.estimate(b, u, p, r) - iter.wsr) / iter.wsr <= iter.tolerance);
					iter.wsr = network.wsr.aggregate.estimate(b, u, p, r);
					iter.counter = iter.counter + 1;
				end
				channel.aggregate.actual = channel_aggregate(channel.direct.actual, channel.forward.actual, channel.backward.actual, reflect.beamformer);
				network.wsr.aggregate.actual(b, u, p, r) = sum(network.weight .* rate_mimo(channel.aggregate.actual, transmit.beamformer, receive.noise), 3);
			end
		end
	end
end
network.wsr.direct.actual = mean(network.wsr.direct.actual, ndims(network.wsr.direct));
network.wsr.aggregate.actual = mean(network.wsr.aggregate.actual, ndims(network.wsr.aggregate.actual));
save('data/ic_wsr_pair.mat');

figure('Name', 'Weighted Sum-Rate vs Number of Transceiver Pairs and Estimation Error', 'Position', [0, 0, 500, 400]);
hold all;
handle.wsr.direct = plot(network.pair, network.wsr.direct.actual / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for u = 1 : number.uncertainty
	handle.wsr.aggregate(1, u) = plot(network.pair, shiftdim(network.wsr.aggregate.actual(1, u, :), 1) / log(2), 'DisplayName', 'D: $\epsilon = ' + string(channel.uncertainty(u)) + '$');
	handle.wsr.aggregate(2, u) = plot(network.pair, shiftdim(network.wsr.aggregate.actual(2, u, :), 1) / log(2), 'DisplayName', 'BD: $\epsilon = ' + string(channel.uncertainty(u)) + '$');
end
style_plot(handle.wsr.aggregate, number.bond);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('Number of Transceiver Pairs');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
savefig('plots/ic_wsr_pair.fig');
matlab2tikz('../assets/simulation/ic_wsr_pair.tex', 'width', '10cm', 'height', '7.5cm');


function [C] = distribution_disk(r, K)
	C = sqrt(r ^ 2 * rand(1, 1, K)) .* exp(1i * 2 * pi * rand(1, 1, K));
	C = [real(C); imag(C)];
end

function [lambda] = pathloss(lambda_0, L, gamma)
	lambda = lambda_0 * L .^ (- gamma);
end

function [W] = precoder_initialize(H, N_e, P_t)
	K = size(H, 3);
	[~, ~, V] = pagesvd(H(:, :, logical(eye(K))));
	W = sqrt(P_t / N_e) * pageswap(V(:, 1 : N_e, :));
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end

function [H] = estimation_effect(H, sigma)
	H = H .* (1 + sigma * sqrt(0.5) * (randn(size(H)) + 1i * randn(size(H))));
end
