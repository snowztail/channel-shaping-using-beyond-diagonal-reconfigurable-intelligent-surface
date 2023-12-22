clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2 .^ (0 : 2 : 4), 256, 4);
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -115));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[transmit.snr, reflect.bond] = deal(pow2db(transmit.power * channel.pathloss.direct ./ receive.noise), [1, reflect.antenna]);
[number.bond, number.noise, number.antenna, number.realization] = deal(length(reflect.bond), length(receive.noise), length(transmit.antenna), 2);

for r = 1 : number.realization
	for a = 1 : number.antenna
		channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna(a));
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna(a));
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
		clear scatter_rate_pc;
		for b = 1 : number.bond
			for n = 1 : number.noise
				transmit.beamformer = precoder_rate_pc(channel.direct, transmit.power, receive.noise(n));
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise(n)));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate_pc(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise(n));
					transmit.beamformer = precoder_rate_pc(channel.aggregate, transmit.power, receive.noise(n));
					receive.rate(n, b, a, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise(n));
					iter.converge = (abs(receive.rate(n, b, a, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate(n, b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.rate = mean(receive.rate, ndims(receive.rate));
save('data/pc_rate_tx.mat');

figure('Name', 'Achievable Rate vs Transmit Antenna', 'Position', [0, 0, 500, 400]);
hold all;
for a = 1 : number.antenna
	for b = 1 : number.bond
		handle.rate(b, a) = plot(transmit.snr, receive.rate(:, b, a) / log(2), 'DisplayName', '$(N^\mathrm{T}, L) = (' + string(transmit.antenna(a)) + ', ' + string(reflect.bond(b)) + ')$');
	end
end
style_plot(handle.rate, number.bond);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Direct SNR [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_tx.fig');
matlab2tikz('../assets/simulation/pc_rate_tx.tex', 'width', '10cm', 'height', '7.5cm');
