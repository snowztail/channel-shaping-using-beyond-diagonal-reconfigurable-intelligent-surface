clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna, network.link] = deal(8, 64, 4, 5);
ris.bond = 2 .^ (0 : 1 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -115));
network.coverage = 1e1;
ris.coordinate = 0;
transmit.coordinate = sqrt(network.coverage ^ 2 * rand(network.link, 1)) .* exp(1i * 2 * pi * rand(network.link, 1));
receive.coordinate = sqrt(network.coverage ^ 2 * rand(network.link, 1)) .* exp(1i * 2 * pi * rand(network.link, 1));
channel.pathloss.direct = db2pow(-30) * abs(transmit.coordinate - receive.coordinate') .^ (-2);
channel.pathloss.forward = db2pow(-30) * abs(transmit.coordinate - ris.coordinate) .^ (-2.4);
channel.pathloss.backward = db2pow(-30) * abs(ris.coordinate - receive.coordinate) .^ (-3);
% channel.pathloss.direct = 1;
% channel.pathloss.forward = 0.01;
% channel.pathloss.backward = 0.01;
% [channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
% channel.rank = min(transmit.antenna, receive.antenna);
channel.rank = 3;
[number.bond, number.noise, number.realization] = deal(length(ris.bond), length(receive.noise), 1e1);

for r = 1 : number.realization
	channel.direct = shiftdim(sqrt(channel.pathloss.direct), -2) .* fading_rayleigh(receive.antenna, transmit.antenna, network.link, network.link);
	channel.forward = shiftdim(sqrt(channel.pathloss.forward), -2) .* fading_rayleigh(ris.antenna, transmit.antenna, 1, network.link);
	channel.backward = shiftdim(sqrt(channel.pathloss.backward), -2) .* fading_rayleigh(receive.antenna, ris.antenna, network.link, 1);

	% * No RIS
	[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-8, 0);
	receive.beamformer = repmat(eye(channel.rank, receive.antenna), [1, 1, network.link, 1]);
	transmit.beamformer = repmat(eye(transmit.antenna, channel.rank), [1, 1, 1, network.link]);
    % for k = 1 : network.link
    %     X = randn(channel.rank, receive.antenna) + 1i * randn(channel.rank, receive.antenna);
    %     X = X * (X' * X) ^ (-0.5);
    %     receiver.beamformer(:, :, k) = X;
    %     X = randn(transmit.antenna, channel.rank) + 1i * randn(transmit.antenna, channel.rank);
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
		receive.beamformer = repmat(eye(channel.rank, receive.antenna), [1, 1, network.link, 1]);
		transmit.beamformer = repmat(eye(transmit.antenna, channel.rank), [1, 1, 1, network.link]);
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
