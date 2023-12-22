clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair, transmit.power, network.coverage] = deal(8, 2 .^ (2 : 2 : 8), 4, 3, 5, db2pow(0), 50);
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = pathloss_disk(network.pair, network.coverage, channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward);
[number.bond, number.antenna, number.realization] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = shiftdim(sqrt(channel.pathloss.direct), -2) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	receive.beamformer = combiner_initial_ic(channel.direct, transmit.stream);
	transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
	[iter.converge, iter.tolerance, iter.counter, iter.interference] = deal(false, 1e-4, 0, interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer)));
	while ~iter.converge
		receive.beamformer = combiner_interference_ic(channel.direct, transmit.beamformer);
		transmit.beamformer = precoder_interference_ic(channel.direct, receive.beamformer, transmit.power);
		receive.interference.direct(r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer));
		iter.converge = (abs(receive.interference.direct(r) - iter.interference) / iter.interference <= iter.tolerance);
		iter.interference = receive.interference.direct(r);
		iter.counter = iter.counter + 1;
	end
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = shiftdim(sqrt(channel.pathloss.forward), -2) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
		channel.backward = shiftdim(sqrt(channel.pathloss.backward), -2) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
		clear scatter_interference_ic;
		for b = 1 : number.bond(a)
			reflect.bond = 2 .^ (b - 1);
			[iter.converge, iter.counter] = deal(false, 0);
			while ~iter.converge
				reflect.beamformer = scatter_interference_ic(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer), pagemtimes(channel.forward, transmit.beamformer), pagemtimes(receive.beamformer, channel.backward), reflect.bond);
				channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer);
				receive.beamformer = combiner_interference_ic(channel.aggregate, transmit.beamformer);
				transmit.beamformer = precoder_interference_ic(channel.aggregate, receive.beamformer, transmit.power);
				receive.interference.aggregate(b, a, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
				iter.converge = (abs(receive.interference.aggregate(b, a, r) - iter.interference) / iter.interference <= iter.tolerance);
				iter.interference = receive.interference.aggregate(b, a, r);
				iter.counter = iter.counter + 1;
			end
		end
	end
end
receive.interference.direct = mean(receive.interference.direct, ndims(receive.interference.direct));
receive.interference.aggregate = mean(receive.interference.aggregate, ndims(receive.interference.aggregate));
save('data/ic_interference_sx.mat');

figure('Name', 'Leakage Interference vs RIS Configuration', 'Position', [0, 0, 500, 400]);
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$');
hold all;
handle.interference.direct = refline(0, receive.interference.direct);
set(handle.interference.direct, 'Color', 'k', 'Marker', 'none', 'DisplayName', '$N^\mathrm{S} = 0$');
for a = 1 : number.antenna
	handle.interference.aggregate(a) = plot(1 : number.bond(a), receive.interference.aggregate(1 : number.bond(a), a), 'DisplayName', '$N^\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
end
style_plot(handle.interference.aggregate);
hold off; grid on; box on; legend('Location', 'ne');
xlabel('RIS Group Size');
ylabel('Leakage Interference [W]');
savefig('plots/ic_interference_sx.fig');
matlab2tikz('../assets/simulation/ic_interference_sx.tex', 'width', '8cm', 'height', '6cm');
