clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 128, 4);
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
% transmit.snr.direct = db2pow(-20 : 5 : 20);
[transmit.snr, reflect.bond] = deal(pow2db(transmit.power .* channel.pathloss.direct ./ receive.noise), [1, 4, reflect.antenna]);
[number.bond, number.power, number.realization, flag.direct] = deal(length(reflect.bond), length(transmit.power), 100, true);

for r = 1 : number.realization
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for p = 1 : number.power
		transmit.beamformer.direct = precoder_rate(channel.direct, transmit.power(p), receive.noise);
		receive.rate.direct(p, r) = rate_mimo(channel.direct, transmit.beamformer.direct, receive.noise);
	end
	% * Have RIS
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	clear scatter_power_max scatter_rate;
	for b = 1 : number.bond
		[reflect.beamformer.decouple, channel.aggregate.decouple] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond(b));
		for p = 1 : number.power
			% * alternating optimization
			transmit.beamformer.alternate = precoder_rate(channel.direct, transmit.power(p), receive.noise);
			[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.alternate, receive.noise));
			while ~iter.converge
				[reflect.beamformer.alternate, channel.aggregate.alternate] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer.alternate, reflect.bond(b), receive.noise);
				transmit.beamformer.alternate = precoder_rate(channel.aggregate.alternate, transmit.power(p), receive.noise);
				receive.rate.alternate(p, b, r) = rate_mimo(channel.aggregate.alternate, transmit.beamformer.alternate, receive.noise);
				iter.converge = (abs(receive.rate.alternate(p, b, r) - iter.rate) / iter.rate <= iter.tolerance);
				iter.rate = receive.rate.alternate(p, b, r);
				iter.counter = iter.counter + 1;
			end

			% * decoupled design
			transmit.beamformer.decouple = precoder_rate(channel.aggregate.decouple, transmit.power(p), receive.noise);
			receive.rate.decouple(p, b, r) = rate_mimo(channel.aggregate.decouple, transmit.beamformer.decouple, receive.noise);
		end
	end
end
receive.rate.direct = mean(receive.rate.direct, ndims(receive.rate.direct));
receive.rate.alternate = mean(receive.rate.alternate, ndims(receive.rate.alternate));
receive.rate.decouple = mean(receive.rate.decouple, ndims(receive.rate.decouple));
save('data/rate_beamforming.mat');

figure('Name', 'Achievable Rate vs Beamforming Design', 'Position', [0, 0, 500, 400]);
hold all;
handle.rate.direct = plot(pow2db(transmit.power), receive.rate.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for b = 1 : number.bond
	handle.rate.aggregate(1, b) = plot(pow2db(transmit.power), receive.rate.alternate(:, b) / log(2), 'DisplayName', 'Alternate: $L = ' + string(reflect.bond(b)) + '$');
	handle.rate.aggregate(2, b) = plot(pow2db(transmit.power), receive.rate.decouple(:, b) / log(2), 'DisplayName', 'Decouple: $L = ' + string(reflect.bond(b)) + '$');
end
style_plot(handle.rate.aggregate, 2);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/rate_beamforming.fig');
matlab2tikz('../assets/simulation/rate_beamforming.tex', 'width', '10cm', 'height', '7.5cm');
