clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair] = deal(8, 2 .^ [-inf, 5, 8], 4, 3, 5);
[transmit.power, receive.noise, receive.weight] = deal(db2pow(-20), db2pow(-75 : -10 : -115), ones(network.pair, 1));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.noise, number.antenna, number.realization] = deal(2, length(receive.noise), length(reflect.antenna), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
		channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
		transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
		transmit.beamformer1 = transmit.beamformer;
		for n = 1 : number.noise
			if reflect.antenna(a) == 0
				% transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
				% transmit.beamformer = precoder_rate_ic_bisection(channel.direct, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
				% F1 = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer, receive.noise(n))
				tic
				transmit.beamformer = precoder_rate_ic(channel.direct, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
				F2 = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer, receive.noise(n))
				toc

				tic
				transmit.beamformer1 = precoder_rate_ic_h(channel.direct, transmit.beamformer1, transmit.power, receive.noise(n), receive.weight);
				F3 = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer1, receive.noise(n))
				toc

				receive.rate(1, n, a, r) = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer, receive.noise(n));
			else
				% reflect.bond = [2, reflect.antenna(a)];
				% for b = 1 : number.bond
				% 	reflect.beamformer = eye(reflect.antenna(a));
				% 	channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer);
				% 	transmit.beamformer = precoder_initial_ic(channel.aggregate, transmit.stream, transmit.power);
				% 	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
				% 	iter.rate = receive.weight' * rate_mimo_ic(channel.aggregate, transmit.beamformer, receive.noise(n));
				% 	while ~iter.converge
					% 	[reflect.beamformer, channel.aggregate] = scatter_rate_ic(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.beamformer, reflect.bond(b), receive.noise(n), receive.weight);
					% 	transmit.beamformer = precoder_rate_ic_bisection(channel.aggregate, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
					% 	receive.rate(b, n, a, r) = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer, receive.noise(n));
					% 	iter.converge = (abs(receive.rate(b, n, a, r) - iter.rate) / iter.rate <= iter.tolerance);
					% 	iter.rate = receive.rate(b, n, a, r);
					% 	iter.counter = iter.counter + 1;
				% 	end
				% end
			end
		end
	end
end
receive.rate = mean(receive.rate, 4);
save('data/ic_rate_sx.mat');

figure('Name', 'Leakage Interference vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
set(gca, 'XLim', [0, log2(max(reflect.antenna))], 'XTick', 0 : log2(max(reflect.antenna)), 'XTickLabel', '$2^' + string(vec(0 : log2(max(reflect.antenna)))) + '$');
for a = 1 : number.antenna
	if reflect.antenna(a) == 0
		handle.interference(a) = refline(0, receive.interference(1, a));
	else
		reflect.bond = 2 .^ (0 : log2(reflect.antenna(a))); number.bond = length(reflect.bond);
		handle.interference(a) = plot(log2(reflect.bond), receive.interference(1 : number.bond, a));
	end
end
hold off; grid on;
style_plot(handle.interference); set(handle.interference(find(reflect.antenna == 0)), 'Color', 'k', 'Marker', 'none');
xlabel('RIS Group Size');
ylabel('Leakage Interference [W]');
legend('$N_s = ' + string(vec(reflect.antenna)) + '$', 'Location', 'sw');
savefig('plots/ic_interference_sx.fig');
