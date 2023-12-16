clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair] = deal(8, 2 .^ [-inf, 5, 8], 4, 3, 5);
[transmit.power, receive.noise, receive.weight] = deal(db2pow(-20), db2pow(-75 : -10 : -115), ones(network.pair, 1));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.noise, number.antenna, number.realization] = deal(2, length(receive.noise), length(reflect.antenna), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	for a = 1 : number.antenna
		if reflect.antenna(a) == 0
			transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
			for n = 1 : number.noise
				transmit.beamformer = precoder_rate_ic(channel.direct, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
				receive.rate(n, 1, a, r) = receive.weight' * rate_mimo_ic(channel.direct, transmit.beamformer, receive.noise(n));
			end
		else
			channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
			channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
			reflect.bond = [1, reflect.antenna(a)];
			for b = 1 : number.bond
				reflect.beamformer = eye(reflect.antenna(a));
				channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer);
				transmit.beamformer = precoder_initial_ic(channel.aggregate, transmit.stream, transmit.power);
				for n = 1 : number.noise
					[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-3, 0);
					iter.rate = receive.weight' * rate_mimo_ic(channel.aggregate, transmit.beamformer, receive.noise(n));
					while ~iter.converge
						[reflect.beamformer, channel.aggregate] = scatter_rate_ic(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.beamformer, reflect.bond(b), receive.noise(n), receive.weight);
						transmit.beamformer = precoder_rate_ic(channel.aggregate, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
						receive.rate(n, b, a, r) = receive.weight' * rate_mimo_ic(channel.aggregate, transmit.beamformer, receive.noise(n));
						iter.converge = (abs(receive.rate(n, b, a, r) - iter.rate) / iter.rate <= iter.tolerance);
						iter.rate = receive.rate(n, b, a, r);
						iter.counter = iter.counter + 1;
					end
				end
			end
		end
	end
end
receive.rate = mean(receive.rate, 4);
save('data/ic_rate_sx.mat');

figure('Name', 'Weighted Sum-Rate vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
for a = 1 : number.antenna
	if reflect.antenna(a) == 0
		handle.rate(1, a) = plot(transmit.snr, receive.rate(:, 1, a) / log(2), 'DisplayName', '$N_s = ' + string(reflect.antenna(a)) + '$');
	else
		reflect.bond = [1, reflect.antenna(a)];
		for b = 1 : number.bond
			handle.rate(b, a) = plot(transmit.snr, receive.rate(:, b, a) / log(2), 'DisplayName', '$(N_s, L) = (' + string(reflect.antenna(a)) + ', ' + string(reflect.bond(b)) + ')$');
		end
	end
end
hold off; grid on;
set(handle.rate(1, find(reflect.antenna == 0)), 'Color', 'k', 'Marker', 'none');
style_plot(handle.rate(:, find(reflect.antenna ~= 0)), number.bond);
xlabel('Direct SNR [dB]');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
legend('Location', 'nw');
savefig('plots/ic_rate_sx.fig');
matlab2tikz('../assets/simulation/ic_rate_sx.tex', 'width', '8cm', 'height', '6cm');
