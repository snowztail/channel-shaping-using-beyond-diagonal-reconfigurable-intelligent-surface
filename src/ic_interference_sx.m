clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 2 .^ (2 : 2 : 8), 4);
[network.coverage, network.pair] = deal(20, 5);
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
[number.bond, number.antenna, number.realization, flag.direct] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 1e2, true);

for r = 1 : number.realization
	[network.coordinate.transmit, network.coordinate.reflect, network.coordinate.receive] = deal(distribution_disk(network.coverage, network.pair), [0; 0], distribution_disk(network.coverage, network.pair));
	channel.pathloss.direct = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.transmit), channel.pathloss.exponent.direct);
	channel.pathloss.forward = pageswap(pathloss(channel.pathloss.reference, vecnorm(network.coordinate.reflect - network.coordinate.transmit), channel.pathloss.exponent.forward));
	channel.pathloss.backward = pathloss(channel.pathloss.reference, vecnorm(network.coordinate.receive - network.coordinate.reflect), channel.pathloss.exponent.backward);
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	channel.interference.direct(r) = interference_leakage(channel.direct);
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
		channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
		for b = 1 : number.bond(a)
			clear scatter_interference;
			reflect.bond = 2 .^ (b - 1);
			[reflect.beamformer, channel.aggregate] = scatter_interference(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.interference.aggregate(b, a, r) = interference_leakage(channel.aggregate);
		end
	end
end
channel.interference.direct = mean(channel.interference.direct, ndims(channel.interference.direct));
channel.interference.aggregate = mean(channel.interference.aggregate, ndims(channel.interference.aggregate));
save('data/ic_interference_sx.mat');

figure('Name', 'Leakage Interference vs RIS Configuration', 'Position', [0, 0, 500, 400]);
if flag.direct
    handle.interference.direct = plot(1 : max(number.bond), repmat(channel.interference.direct, [1, max(number.bond)]), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
	hold on;
end
for a = 1 : number.antenna
	handle.interference.aggregate(a) = plot(1 : number.bond(a), channel.interference.aggregate(1 : number.bond(a), a), 'DisplayName', '$N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
end
style_plot(handle.interference.aggregate);
set(handle.interference.aggregate, {'Marker'}, {'none'});
if flag.direct
	legend([handle.interference.direct, handle.interference.aggregate], 'Location', 'ne');
else
	legend('Location', 'ne');
end
hold off; grid on; box on;
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$');
xlabel('RIS Group Size');
ylabel('Leakage Interference [W]');
savefig('plots/ic_interference_sx.fig');
matlab2tikz('../assets/simulation/ic_interference_sx.tex', 'width', '8cm', 'height', '6cm');


function [C] = distribution_disk(r, K)
	C = sqrt(r ^ 2 * rand(1, 1, K)) .* exp(1i * 2 * pi * rand(1, 1, K));
	C = [real(C); imag(C)];
end

function [lambda] = pathloss(lambda_0, L, gamma)
	lambda = lambda_0 * L .^ (- gamma);
end

function [A] = pageswap(A)
	A = permute(A, [1, 2, 4, 3]);
end
