clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair] = deal(8, 2 .^ [5, 8], 4, 3, 5);
[transmit.power, receive.noise, receive.weight] = deal(db2pow(-20), db2pow(-75 : -10 : -115), ones(network.pair, 1));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.noise, number.antenna, number.realization] = deal(3, length(receive.noise), length(reflect.antenna), 2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
	for n = 1 : number.noise
		transmit.beamformer = precoder_rate_explicit_ic(channel.direct, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
		receive.wsr.direct(n, r) = receive.weight' * rate_mimo(channel.direct, transmit.beamformer, receive.noise(n));
	end
	% * Have RIS
	for a = 1 : number.antenna
		reflect.bond = [1, 4, reflect.antenna(a)];
		channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
		channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
		clear scatter_rate_ic;
		for b = 1 : number.bond
			transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
			for n = 1 : number.noise
				[iter.converge, iter.tolerance, iter.counter, iter.wsr] = deal(false, 1e-4, 0, receive.weight' * rate_mimo(channel.direct, transmit.beamformer, receive.noise(n)));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate_ic(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise(n), receive.weight);
					transmit.beamformer = precoder_rate_explicit_ic(channel.aggregate, transmit.beamformer, transmit.power, receive.noise(n), receive.weight);
					receive.wsr.aggregate(n, b, a, r) = receive.weight' * rate_mimo(channel.aggregate, transmit.beamformer, receive.noise(n));
					iter.converge = (abs(receive.wsr.aggregate(n, b, a, r) - iter.wsr) / iter.wsr <= iter.tolerance);
					iter.wsr = receive.wsr.aggregate(n, b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.wsr.direct = mean(receive.wsr.direct, ndims(receive.wsr.direct));
receive.wsr.aggregate = mean(receive.wsr.aggregate, ndims(receive.wsr.aggregate));
save('data/ic_rate_sx.mat');

figure('Name', 'Weighted Sum-Rate vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
handle.wsr.direct = plot(transmit.snr, receive.wsr.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', '$N^\mathrm{S} = 0$');
for a = 1 : number.antenna
	reflect.bond = [1, 4, reflect.antenna(a)];
	for b = 1 : number.bond
		handle.wsr.aggregate(b, a) = plot(transmit.snr, receive.wsr.aggregate(:, b, a) / log(2), 'DisplayName', '$(N^\mathrm{S}, L) = (' + string(reflect.antenna(a)) + ', ' + string(reflect.bond(b)) + ')$');
	end
end
style_plot(handle.wsr.aggregate, number.bond);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('Direct SNR [dB]');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
savefig('plots/ic_rate_sx.fig');
matlab2tikz('../assets/simulation/ic_rate_sx.tex', 'width', '10cm', 'height', '7.5cm');
