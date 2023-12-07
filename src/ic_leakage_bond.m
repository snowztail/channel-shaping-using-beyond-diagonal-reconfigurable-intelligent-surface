clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna, transmit.stream, network.pair, network.coverage] = deal(8, 2 .^ [-inf, 0 : 2 : 8], 4, 3, 5, 50);
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -115));
[channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(db2pow(-30), 3, 2.4, 2.4);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = pathloss_disk(network.pair, network.coverage, channel.pathloss.reference, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward);
[number.antenna, number.noise, number.realization] = deal(length(ris.antenna), length(receive.noise), 1e1);

for r = 1 : number.realization
	channel.direct = shiftdim(sqrt(channel.pathloss.direct), -2) .* fading_nlos(receive.antenna, transmit.antenna, network.pair, network.pair);
	for a = 1 : number.antenna
		ris.bond = 2 .^ [-inf, (0 : log2(ris.antenna(a)))]; number.bond = length(ris.bond);
		channel.forward = shiftdim(sqrt(channel.pathloss.forward), -2) .* fading_nlos(ris.antenna(a), transmit.antenna, 1, network.pair);
		channel.backward = shiftdim(sqrt(channel.pathloss.backward), -2) .* fading_nlos(receive.antenna, ris.antenna(a), network.pair, 1);
		for b = 1 : number.bond
			[iter.converge, iter.tolerance, iter.counter, iter.leakage] = deal(false, 1e-8, 0, 0);
			ris.scatter = eye(ris.antenna(a));
			channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
			transmit.beamformer = precoder_initialize_ic(channel.aggregate, transmit.stream);
			receive.beamformer = combiner_leakage_ic(channel.aggregate, transmit.beamformer);
			receive.leakage(b, a, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
			while ~iter.converge
				iter.leakage = receive.leakage(b, a, r);
				receive.beamformer = combiner_leakage_ic(channel.aggregate, transmit.beamformer);
				transmit.beamformer = precoder_leakage_ic(channel.aggregate, receive.beamformer);
				ris.scatter = scatter_leakage_ic(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer), pagemtimes(channel.forward, transmit.beamformer), pagemtimes(receive.beamformer, channel.backward), ris.scatter, ris.bond(b));
				channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
				receive.leakage(b, a, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
				iter.converge = (abs(receive.leakage(b, a, r) - iter.leakage) <= iter.tolerance);
				iter.counter = iter.counter + 1;
			end
		end
	end



	ris.bond = 2 .^ (0 : 1 : log2(ris.antenna));

	% * No RIS
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-8, 0);


	receive.beamformer = repmat(eye(transmit.stream, receive.antenna), [1, 1, network.pair, 1]);
	transmit.beamformer = repmat(eye(transmit.antenna, transmit.stream), [1, 1, 1, network.pair]);
	% for k = 1 : network.pair
	%     X = randn(transmit.stream, receive.antenna) + 1i * randn(transmit.stream, receive.antenna);
	%     X = X * (X' * X) ^ (-0.5);
	%     receiver.beamformer(:, :, k) = X;
	%     X = randn(transmit.antenna, transmit.stream) + 1i * randn(transmit.antenna, transmit.stream);
	%     X = X * (X' * X) ^ (-0.5);
	%     transmit.beamformer(:, :, :, k) = X;
	% end
	receive.interference.direct(r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer));
	while ~iter.converge
		iter.interference = receive.interference.direct(r);
		receive.beamformer = receive_min_interference(channel.direct, transmit.beamformer);
		transmit.beamformer = transmit_min_interference(channel.direct, receive.beamformer);
		receive.interference.direct(r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.direct), transmit.beamformer));
		iter.converge = (abs(receive.interference.direct(r) - iter.interference) <= iter.tolerance);
		iter.counter = iter.counter + 1;
	end
	for n = 1 : number.noise
		receive.rate.direct(n, r) = sum(rate_ic(channel.direct, transmit.beamformer, receive.beamformer, receive.noise(n)), 'all');
	end

	% * With RIS
	for b = 1 : number.bond
		[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-8, 0);
		receive.beamformer = repmat(eye(transmit.stream, receive.antenna), [1, 1, network.pair, 1]);
		transmit.beamformer = repmat(eye(transmit.antenna, transmit.stream), [1, 1, 1, network.pair]);
		ris.scatter = eye(ris.antenna);
		channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
		receive.interference.aggregate(b, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
		while ~iter.converge
			iter.interference = receive.interference.aggregate(b, r);
			receive.beamformer = receive_min_interference(channel.aggregate, transmit.beamformer);
			transmit.beamformer = transmit_min_interference(channel.aggregate, receive.beamformer);
			[channel.equivalent.direct, channel.equivalent.forward, channel.equivalent.backward] = channel_equivalent(channel.direct, channel.forward, channel.backward, transmit.beamformer, receive.beamformer);
			ris.scatter = ris_min_interference(channel.equivalent.direct, channel.equivalent.forward, channel.equivalent.backward, ris.scatter, ris.group(b));
			channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
			receive.interference.aggregate(b, r) = interference_leakage(pagemtimes(pagemtimes(receive.beamformer, channel.aggregate), transmit.beamformer));
			iter.converge = (abs(receive.interference.aggregate(b, r) - iter.interference) <= iter.tolerance);
			iter.counter = iter.counter + 1;
		end
		for n = 1 : number.noise
			receive.rate.aggregate(b, n, r) = sum(rate_ic(channel.aggregate, transmit.beamformer, receive.beamformer, receive.noise(n)), 'all');
		end
	end
end

receive.interference.direct = mean(receive.interference.direct, 2);
receive.interference.aggregate = mean(receive.interference.aggregate, 2);
receive.rate.direct = mean(receive.rate.direct, 2);
receive.rate.aggregate = mean(receive.rate.aggregate, 3);

figure('Name', 'Total Leakage Interference vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.interference(1) = scatter(0, receive.interference.direct, 'Marker', 'o', 'DisplayName', 'No RIS');
hold on;
handle.interference(2) = plot(ris.bond, receive.interference.aggregate, 'Marker', 'x', 'DisplayName', 'BD-RIS');
legend(handle.interference, 'Location', 'nw'); grid on; box on; axis tight;
xlabel('RIS Group Size');
ylabel('Total Interference Leakage [W]');
savefig('plots/min_interference_bound_il.fig');
% savefig('plots/min_interference_bound_il_circle.fig');

figure('Name', 'Total Rate vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.rate = gobjects(number.bond + 1, 1);
hold all;
handle.rate(1) = plot(pow2db(receive.noise), receive.rate.direct / log(2), 'DisplayName', 'No RIS');
for b = 1 : number.bond
	handle.rate(b + 1) = plot(pow2db(receive.noise), receive.rate.aggregate(b, :) / log(2), 'DisplayName', strcat('$N_g = ', num2str(ris.bond(b)), '$'));
end
hold off; legend('Location', 'nw'); grid on; box on; axis tight;
xlabel('Average Noise Power [dB]');
ylabel('Total Rate [bits/s/Hz]');
style_plot(handle.rate);
% savefig('plots/min_interference_bound_rate.fig');
savefig('plots/min_interference_bound_rate_circle.fig');
