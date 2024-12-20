clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 128, 4);
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[channel.kfactor.direct, channel.kfactor.forward, channel.kfactor.backward] = deal(0, [0, 10, inf], [0, 10, inf]);
[transmit.snr, reflect.bond] = deal(pow2db(transmit.power * channel.pathloss.direct ./ receive.noise), [1, reflect.antenna]);
[number.bond, number.power, number.kfactor, number.realization] = deal(length(reflect.bond), length(transmit.power), length(channel.kfactor.forward), 1e2);

for r = 1 : number.realization
	for k = 1 : number.kfactor
		channel.direct = sqrt(channel.pathloss.direct) * fading_ricean(receive.antenna, transmit.antenna, channel.kfactor.direct);
		channel.forward = sqrt(channel.pathloss.forward) * fading_ricean(reflect.antenna, transmit.antenna, channel.kfactor.forward(k));
		channel.backward = sqrt(channel.pathloss.backward) * fading_ricean(receive.antenna, reflect.antenna, channel.kfactor.backward(k));
		clear scatter_rate;
		for b = 1 : number.bond
			for p = 1 : number.power
				transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise);
					transmit.beamformer = precoder_rate(channel.aggregate, transmit.power(p), receive.noise);
					receive.rate(p, b, k, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise);
					iter.converge = (abs(receive.rate(p, b, k, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate(p, b, k, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.rate = mean(receive.rate, ndims(receive.rate));
save('data/pc_rate_kfactor.mat');

figure('Name', 'Achievable Rate vs Ricean K Factor', 'Position', [0, 0, 500, 400]);
hold all;
for k = 1 : number.kfactor
	if isinf(channel.kfactor.forward(k))
		handle.rate(1, k) = plot(pow2db(transmit.power), receive.rate(:, 1, k) / log(2), 'DisplayName', '{D: $\kappa_\mathrm{F} = \kappa_\mathrm{B} \to \infty$}');
		handle.rate(2, k) = plot(pow2db(transmit.power), receive.rate(:, 2, k) / log(2), 'DisplayName', '{BD: $\kappa_\mathrm{F} = \kappa_\mathrm{B} \to \infty$}');
	else
		handle.rate(1, k) = plot(pow2db(transmit.power), receive.rate(:, 1, k) / log(2), 'DisplayName', 'D: $\kappa_\mathrm{F} = \kappa_\mathrm{B} = ' + string(channel.kfactor.forward(k)) + '$');
		handle.rate(2, k) = plot(pow2db(transmit.power), receive.rate(:, 2, k) / log(2), 'DisplayName', 'BD: $\kappa_\mathrm{F} = \kappa_\mathrm{B} = ' + string(channel.kfactor.forward(k)) + '$');
	end
end
style_plot(handle.rate, 2);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_kfactor.fig');
matlab2tikz('../assets/simulation/pc_rate_kfactor.tex', 'width', '10cm', 'height', '7.5cm');
