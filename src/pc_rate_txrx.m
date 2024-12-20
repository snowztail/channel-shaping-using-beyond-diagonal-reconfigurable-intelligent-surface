clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2 .^ (0 : 2 : 4), 128, 2 .^ (0 : 2 : 4));
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[transmit.snr, reflect.bond] = deal(pow2db(transmit.power * channel.pathloss.direct ./ receive.noise), [1, reflect.antenna]);
[number.bond, number.power, number.antenna, number.realization] = deal(length(reflect.bond), length(transmit.power), length(transmit.antenna), 1e2);

for r = 1 : number.realization
	for a = 1 : number.antenna
		channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna(a), transmit.antenna(a));
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna(a));
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna(a), reflect.antenna);
		clear scatter_rate;
		for b = 1 : number.bond
			for p = 1 : number.power
				transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise);
					transmit.beamformer = precoder_rate(channel.aggregate, transmit.power(p), receive.noise);
					receive.rate(p, b, a, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise);
					iter.converge = (abs(receive.rate(p, b, a, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate(p, b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.rate = mean(receive.rate, ndims(receive.rate));
save('data/pc_rate_txrx.mat');

figure('Name', 'Achievable Rate vs Transmit and Receive Antenna', 'Position', [0, 0, 500, 400]);
hold all;
for a = 1 : number.antenna
	handle.rate(1, a) = plot(pow2db(transmit.power), receive.rate(:, 1, a) / log(2), 'DisplayName', 'D: $N_\mathrm{T}=N_\mathrm{R} = ' + string(transmit.antenna(a)) + '$');
	handle.rate(2, a) = plot(pow2db(transmit.power), receive.rate(:, 2, a) / log(2), 'DisplayName', 'BD: $N_\mathrm{T}=N_\mathrm{R} = ' + string(transmit.antenna(a)) + '$');
end
style_plot(handle.rate, 2);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_txrx.fig');
matlab2tikz('../assets/simulation/pc_rate_txrx.tex', 'width', '10cm', 'height', '7.5cm');
