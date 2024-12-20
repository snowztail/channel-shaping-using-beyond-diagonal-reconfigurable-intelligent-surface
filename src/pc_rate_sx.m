clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 2 .^ (4 : 2 : 8), 4);
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.power, number.antenna, number.realization] = deal(2, length(transmit.power), length(reflect.antenna), 2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for p = 1 : number.power
		transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
		receive.rate.direct(p, r) = rate_mimo(channel.direct, transmit.beamformer, receive.noise);
	end
	% * Have RIS
	for a = 1 : number.antenna
		reflect.bond = [1, reflect.antenna(a)];
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		clear scatter_rate;
		for b = 1 : number.bond
			for p = 1 : number.power
				transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise);
					transmit.beamformer = precoder_rate(channel.aggregate, transmit.power(p), receive.noise);
					receive.rate.aggregate(p, b, a, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise);
					iter.converge = (abs(receive.rate.aggregate(p, b, a, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate.aggregate(p, b, a, r);
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
handle.rate.direct = plot(pow2db(transmit.power), receive.rate.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for a = 1 : number.antenna
	handle.rate.aggregate(1, a) = plot(pow2db(transmit.power), receive.rate.aggregate(:, 1, a) / log(2), 'DisplayName', 'D: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
	handle.rate.aggregate(2, a) = plot(pow2db(transmit.power), receive.rate.aggregate(:, 2, a) / log(2), 'DisplayName', 'BD: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
end
style_plot(handle.rate.aggregate, number.bond);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_sx.fig');
matlab2tikz('../assets/simulation/pc_rate_sx.tex', 'width', '10cm', 'height', '7.5cm');
