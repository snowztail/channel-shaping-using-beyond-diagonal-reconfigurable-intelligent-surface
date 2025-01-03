clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream] = deal(2, 128, 2, 2);
[network.coverage, network.pair] = deal(20, 2);
[transmit.power, receive.noise, network.weight] = deal(db2pow(-20 : 5 : 20), db2pow(-75), ones(1, 1, network.pair));
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
reflect.bond = [1, reflect.antenna];
[number.bond, number.power, number.realization] = deal(2, length(transmit.power), 1e2);

for r = 1 : number.realization
	[network.coordinate.transmit, network.coordinate.reflect, network.coordinate.receive] = deal(distribution_disk(network.coverage, network.pair), [0; 0], distribution_disk(network.coverage, network.pair));
	channel.pathloss.direct = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.transmit), channel.pathloss.exponent.direct);
	channel.pathloss.forward = pageswap(pathloss(channel.pathloss.reference, vecnorm(network.coordinate.reflect - network.coordinate.transmit), channel.pathloss.exponent.forward));
	channel.pathloss.backward = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.reflect), channel.pathloss.exponent.backward);
	channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna, transmit.antenna, 1, network.pair);
	channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna, network.pair, 1);
	% * No RIS
	for p = 1 : number.power
		transmit.beamformer.direct = precoder_initialize(channel.direct, transmit.stream, transmit.power(p));
		transmit.beamformer.direct = precoder_wsr(channel.direct, transmit.beamformer.direct, transmit.power(p), receive.noise, network.weight);
		network.wsr.direct(p, r) = sum(network.weight .* rate_mimo(channel.direct, transmit.beamformer.direct, receive.noise), 3);
	end
	% * Have RIS
    clear scatter_interference;
	for b = 1 : number.bond
		[reflect.beamformer.decouple, channel.aggregate.decouple] = scatter_interference(channel.direct, channel.forward, channel.backward, reflect.bond(b));
		clear scatter_wsr;
		for p = 1 : number.power
			% * alternating optimization
			transmit.beamformer.alternate = precoder_initialize(channel.direct, transmit.stream, transmit.power(p));
			% transmit.beamformer.alternate = precoder_wsr(channel.direct, transmit.beamformer.alternate, transmit.power(p), receive.noise, network.weight);
			[iter.converge, iter.tolerance, iter.counter, iter.wsr] = deal(false, 1e-3, 0, sum(network.weight .* rate_mimo(channel.direct, transmit.beamformer.alternate, receive.noise), 3));
			while ~iter.converge && iter.counter <= 1e2
				[reflect.beamformer.alternate, channel.aggregate.alternate] = scatter_wsr(channel.direct, channel.forward, channel.backward, transmit.beamformer.alternate, reflect.bond(b), receive.noise, network.weight);
				transmit.beamformer.alternate = precoder_wsr(channel.aggregate.alternate, transmit.beamformer.alternate, transmit.power(p), receive.noise, network.weight);
				network.wsr.alternate(p, b, r) = sum(network.weight .* rate_mimo(channel.aggregate.alternate, transmit.beamformer.alternate, receive.noise), 3);
				iter.converge = (abs(network.wsr.alternate(p, b, r) - iter.wsr) / iter.wsr <= iter.tolerance);
				iter.wsr = network.wsr.alternate(p, b, r);
				iter.counter = iter.counter + 1;
			end

			% * decoupled design
			transmit.beamformer.decouple = precoder_initialize(channel.direct, transmit.stream, transmit.power(p));
			transmit.beamformer.decouple = precoder_wsr(channel.aggregate.decouple, transmit.beamformer.decouple, transmit.power(p), receive.noise, network.weight);
			network.wsr.decouple(p, b, r) = sum(network.weight .* rate_mimo(channel.aggregate.decouple, transmit.beamformer.decouple, receive.noise), 3);
		end
	end
end
network.wsr.direct = mean(network.wsr.direct, ndims(network.wsr.direct));
network.wsr.alternate = mean(network.wsr.alternate, ndims(network.wsr.alternate));
network.wsr.decouple = mean(network.wsr.decouple, ndims(network.wsr.decouple));
save('data/ic_wsr_beamforming.mat');

figure('Name', 'Weighted Sum-Rate vs Beamforming Design', 'Position', [0, 0, 500, 400]);
hold all;
handle.wsr.direct = plot(pow2db(transmit.power), network.wsr.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for b = 1 : number.bond
	handle.wsr.aggregate(1, b) = plot(pow2db(transmit.power), network.wsr.alternate(:, b) / log(2), 'DisplayName', 'Alternate: $L = ' + string(reflect.bond(b)) + '$');
	handle.wsr.aggregate(2, b) = plot(pow2db(transmit.power), network.wsr.decouple(:, b) / log(2), 'DisplayName', 'Decouple: $L = ' + string(reflect.bond(b)) + '$');
end
style_plot(handle.wsr.aggregate, number.bond);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
savefig('plots/ic_wsr_beamforming.fig');
matlab2tikz('../assets/simulation/ic_wsr_beamforming.tex', 'width', '10cm', 'height', '7.5cm');


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
