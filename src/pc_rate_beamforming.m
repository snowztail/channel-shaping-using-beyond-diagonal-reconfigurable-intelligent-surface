clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2, 64, 2);
[transmit.power, receive.noise] = deal(db2pow(-20 : 10 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[transmit.snr, reflect.bond] = deal(pow2db(transmit.power .* channel.pathloss.direct ./ receive.noise), [1, reflect.antenna]);
[number.bond, number.power, number.realization, flag.direct] = deal(length(reflect.bond), length(transmit.power), 2, true);

for r = 1 : number.realization
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	clear scatter_power_max_pc scatter_rate_pc;
	for b = 1 : number.bond
		[reflect.beamformer.decouple, channel.aggregate.decouple] = scatter_power_max_pc(channel.direct, channel.forward, channel.backward, reflect.bond(b));
		for p = 1 : number.power
			% * decoupled design
			transmit.beamformer.decouple = precoder_rate_pc(channel.aggregate.decouple, transmit.power(p), receive.noise);
			receive.rate.decouple(p, b, r) = rate_mimo(channel.aggregate.decouple, transmit.beamformer.decouple, receive.noise);

			% * alternating optimization
			transmit.beamformer.alternate = precoder_rate_pc(channel.direct, transmit.power(p), receive.noise);
			[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.alternate, receive.noise));
			while ~iter.converge
				[reflect.beamformer.alternate, channel.aggregate.alternate] = scatter_rate_pc(channel.direct, channel.forward, channel.backward, transmit.beamformer.alternate, reflect.bond(b), receive.noise);
				transmit.beamformer.alternate = precoder_rate_pc(channel.aggregate.alternate, transmit.power(p), receive.noise);
				receive.rate.alternate(p, b, r) = rate_mimo(channel.aggregate.alternate, transmit.beamformer.alternate, receive.noise);
				iter.converge = (abs(receive.rate.alternate(p, b, r) - iter.rate) / iter.rate <= iter.tolerance);
				iter.rate = receive.rate.alternate(p, b, r);
				iter.counter = iter.counter + 1;
			end
		end
	end
end
receive.rate.decouple = mean(receive.rate.decouple, ndims(receive.rate.decouple));
receive.rate.alternate = mean(receive.rate.alternate, ndims(receive.rate.alternate));
save('data/pc_rate_beamforming.mat');

figure('Name', 'Achievable Rate vs Beamforming Design', 'Position', [0, 0, 500, 400]);
hold all;
for b = 1 : number.bond
	handle.rate.decouple(b) = plot(pow2db(transmit.power), receive.rate.decouple(:, b) / log(2), 'DisplayName', 'Decouple: $L = ' + string(reflect.bond(b)) + '$');
	handle.rate.alternate(b) = plot(pow2db(transmit.power), receive.rate.alternate(:, b) / log(2), 'DisplayName', 'Alternate: $L = ' + string(reflect.bond(b)) + '$');
end
style_plot(handle.rate.decouple, number.bond);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_beamforming.fig');
matlab2tikz('../assets/simulation/pc_rate_beamforming.tex', 'width', '10cm', 'height', '7.5cm');
