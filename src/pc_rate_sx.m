clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(8, 2 .^ [5, 8], 4);
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -115));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.noise, number.antenna, number.realization] = deal(3, length(receive.noise), length(reflect.antenna), 2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for n = 1 : number.noise
		transmit.beamformer = precoder_rate_pc(channel.direct, transmit.power, receive.noise(n));
		receive.rate.direct(n, r) = rate_mimo(channel.direct, transmit.beamformer, receive.noise(n));
	end
	% * Have RIS
	for a = 1 : number.antenna
		reflect.bond = [1, 4, reflect.antenna(a)];
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		clear scatter_rate_pc;
		for b = 1 : number.bond
			for n = 1 : number.noise
				transmit.beamformer = precoder_rate_pc(channel.direct, transmit.power, receive.noise(n));
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise(n)));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate_pc(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise(n));
					transmit.beamformer = precoder_rate_pc(channel.aggregate, transmit.power, receive.noise(n));
					receive.rate.aggregate(n, b, a, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise(n));
					iter.converge = (abs(receive.rate.aggregate(n, b, a, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate.aggregate(n, b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.rate.direct = mean(receive.rate.direct, ndims(receive.rate.direct));
receive.rate.aggregate = mean(receive.rate.aggregate, ndims(receive.rate.aggregate));
save('data/pc_rate_sx.mat');

figure('Name', 'Achievable Rate vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
handle.rate.direct = plot(transmit.snr, receive.rate.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', '$N^\mathrm{S} = 0$');
for a = 1 : number.antenna
	reflect.bond = [1, 4, reflect.antenna(a)];
	for b = 1 : number.bond
		handle.rate.aggregate(b, a) = plot(transmit.snr, receive.rate.aggregate(:, b, a) / log(2), 'DisplayName', '$(N^\mathrm{S}, L) = (' + string(reflect.antenna(a)) + ', ' + string(reflect.bond(b)) + ')$');
	end
end
style_plot(handle.rate.aggregate, number.bond);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Direct SNR [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_sx.fig');
matlab2tikz('../assets/simulation/pc_rate_sx.tex', 'width', '10cm', 'height', '7.5cm');
