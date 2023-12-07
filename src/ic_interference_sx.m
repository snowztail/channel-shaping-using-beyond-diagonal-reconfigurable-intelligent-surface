clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair, transmit.power, network.coverage] = deal(8, 2 .^ [-inf, 2 : 2 : 8], 4, 3, 5, db2pow(0), 50);
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = pathloss_disk(network.pair, network.coverage, channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward);
[number.antenna, number.realization] = deal(length(reflect.antenna), 1e1);

for r = 1 : number.realization
	channel.direct = shiftdim(sqrt(channel.pathloss.direct), -2) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	for a = 1 : number.antenna
		if reflect.antenna(a) == 0
			receive.beamformer = combiner_initial_ic(channel.direct, transmit.stream);
			transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
			[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
			iter.interference = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer));
			while ~iter.converge
				receive.beamformer = combiner_interference_ic(channel.direct, transmit.beamformer);
				transmit.beamformer = precoder_interference_ic(channel.direct, receive.beamformer, transmit.power);
				receive.interference(1, a, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer));
				iter.converge = (abs(receive.interference(1, a, r) - iter.interference) / iter.interference <= iter.tolerance);
				iter.interference = receive.interference(1, a, r);
				iter.counter = iter.counter + 1;
			end
		else
			channel.forward = shiftdim(sqrt(channel.pathloss.forward), -2) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair);
			channel.backward = shiftdim(sqrt(channel.pathloss.backward), -2) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair, 1);
			reflect.bond = 2 .^ (0 : log2(reflect.antenna(a))); number.bond = length(reflect.bond);
			for b = 1 : number.bond
				reflect.beamformer = eye(reflect.antenna(a));
				channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer);
				receive.beamformer = combiner_initial_ic(channel.aggregate, transmit.stream);
				transmit.beamformer = precoder_initial_ic(channel.aggregate, transmit.stream, transmit.power);
				[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
				iter.interference = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
				while ~iter.converge
					reflect.beamformer = scatter_interference_ic(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer), pagemtimes(channel.forward, transmit.beamformer), pagemtimes(receive.beamformer, channel.backward), reflect.beamformer, reflect.bond(b));
					channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer);
					receive.beamformer = combiner_interference_ic(channel.aggregate, transmit.beamformer);
					transmit.beamformer = precoder_interference_ic(channel.aggregate, receive.beamformer, transmit.power);
					receive.interference(b, a, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
					iter.converge = (abs(receive.interference(b, a, r) - iter.interference) / iter.interference <= iter.tolerance);
					iter.interference = receive.interference(b, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.interference = mean(receive.interference, 3);
save('data/ic_interference_sx.mat');

figure('Name', 'Leakage Interference vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
set(gca, 'XLim', [0, log2(max(reflect.antenna))], 'XTick', 0 : log2(max(reflect.antenna)), 'XTickLabel', '$2^' + string(vec(0 : log2(max(reflect.antenna)))) + '$');
for a = 1 : number.antenna
	if reflect.antenna(a) == 0
		handle.interference(a) = refline(0, receive.interference(1, a));
	else
		reflect.bond = 2 .^ (0 : log2(reflect.antenna(a))); number.bond = length(reflect.bond);
		handle.interference(a) = plot(log2(reflect.bond), receive.interference(1 : number.bond, a));
	end
end
hold off; grid on;
style_plot(handle.interference); set(handle.interference(find(reflect.antenna == 0)), 'Color', 'k', 'Marker', 'none');
xlabel('RIS Group Size');
ylabel('Leakage Interference [W]');
legend('$N_s = ' + string(vec(reflect.antenna)) + '$', 'Location', 'sw');
savefig('plots/ic_interference_sx.fig');
