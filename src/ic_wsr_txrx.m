clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream] = deal(2 .^ (2 : 2 : 4), 128, 2 .^ (2 : 2 : 4), 2);
[network.coverage, network.pair] = deal(20, 10);
[transmit.power, receive.noise, network.weight] = deal(db2pow(-20 : 5 : 20), db2pow(-75), ones(1, 1, network.pair));
[network.coordinate.transmit, network.coordinate.reflect, network.coordinate.receive] = deal(distribution_disk(network.coverage, network.pair), [0; 0], distribution_disk(network.coverage, network.pair));
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
channel.pathloss.direct = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.transmit), channel.pathloss.exponent.direct);
channel.pathloss.forward = pageswap(pathloss(channel.pathloss.reference, vecnorm(network.coordinate.reflect - network.coordinate.transmit), channel.pathloss.exponent.forward));
channel.pathloss.backward = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.reflect), channel.pathloss.exponent.backward);
reflect.bond = [1, reflect.antenna];
[number.bond, number.power, number.antenna, number.realization] = deal(2, length(transmit.power), length(transmit.antenna), 1e2);

for r = 1 : number.realization
	for a = 1 : number.antenna
		channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna(a), transmit.antenna(a), network.pair, network.pair);
		channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna, transmit.antenna(a), 1, network.pair);
		channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna(a), reflect.antenna, network.pair, 1);
		clear scatter_wsr;
		for b = 1 : number.bond
			for p = 1 : number.power
				transmit.beamformer = precoder_initialize(channel.direct, transmit.stream, transmit.power(p));
				[iter.converge, iter.tolerance, iter.counter, iter.wsr] = deal(false, 1e-4, 0, sum(network.weight .* rate_mimo(channel.direct, transmit.beamformer, receive.noise), 3));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_wsr(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise, network.weight);
					transmit.beamformer = precoder_wsr(channel.aggregate, transmit.beamformer, transmit.power(p), receive.noise, network.weight);
					network.wsr.aggregate(p, b, a, r) = sum(network.weight .* rate_mimo(channel.aggregate, transmit.beamformer, receive.noise), 3);
					iter.converge = (abs(network.wsr.aggregate(p, b, a, r) - iter.wsr) / iter.wsr <= iter.tolerance);
					iter.wsr = network.wsr.aggregate(p, b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
network.wsr.aggregate = mean(network.wsr.aggregate, ndims(network.wsr.aggregate));
save('data/ic_wsr_txrx.mat');

figure('Name', 'Weighted Sum-Rate vs Transmit and Receive Antenna', 'Position', [0, 0, 500, 400]);
hold all;
for a = 1 : number.antenna
	handle.wsr.aggregate(1, a) = plot(pow2db(transmit.power), network.wsr.aggregate(:, 1, a) / log(2), 'DisplayName', 'D: $N_\mathrm{T}=N_\mathrm{R} = ' + string(transmit.antenna(a)) + '$');
	handle.wsr.aggregate(2, a) = plot(pow2db(transmit.power), network.wsr.aggregate(:, 2, a) / log(2), 'DisplayName', 'BD: $N_\mathrm{T}=N_\mathrm{R} = ' + string(transmit.antenna(a)) + '$');
end
style_plot(handle.wsr.aggregate, number.bond);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
savefig('plots/ic_wsr_txrx.fig');
matlab2tikz('../assets/simulation/ic_wsr_txrx.tex', 'width', '10cm', 'height', '7.5cm');


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
