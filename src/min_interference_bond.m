clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna, network.link] = deal(2, 256, 4, 3);
ris.bond = 2 .^ (0 : 2 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -105));
network.coverage = 1e2;
ris.coordinate = 0;
transmit.coordinate = sqrt(network.coverage ^ 2 * rand(network.link, 1)) .* exp(1i * 2 * pi * rand(network.link, 1));
receiver.coordinate = sqrt(network.coverage ^ 2 * rand(network.link, 1)) .* exp(1i * 2 * pi * rand(network.link, 1));
channel.pathloss.direct = db2pow(-30) * abs(transmit.coordinate - receiver.coordinate') .^ (-2);
channel.pathloss.forward = db2pow(-30) * abs(transmit.coordinate - ris.coordinate) .^ (-2.4);
channel.pathloss.backward = db2pow(-30) * abs(ris.coordinate - receiver.coordinate) .^ (-3);
channel.rank = min(transmit.antenna, receive.antenna);
[number.bond, number.realization] = deal(length(ris.bond), 1e1);

for r = 1 : number.realization
	channel.direct = shiftdim(sqrt(channel.pathloss.direct), -2) .* fading_rayleigh(receive.antenna, transmit.antenna, network.link, network.link);
	channel.forward = shiftdim(sqrt(channel.pathloss.forward), -2) .* fading_rayleigh(ris.antenna, transmit.antenna, network.link);
	channel.backward = shiftdim(sqrt(channel.pathloss.backward), -2) .* fading_rayleigh(receive.antenna, ris.antenna, network.link);
	for b = 1 : number.bond
		[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-5, 0);
		receive.beamformer = repmat(eye(channel.rank, receive.antenna), [1, 1, network.link]);
		transmit.beamformer = repmat(eye(transmit.antenna, channel.rank), [1, 1, network.link]);
		ris.scatter = eye(ris.antenna);
		channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
		
		channel.equivalent.direct = pagemtimes(permute(receive.beamformer, [1, 2, 4, 3]), pagemtimes(channel.direct, transmit.beamformer));
		channel.equivalent.forward = pagemtimes(channel.forward, transmit.beamformer);
		channel.equivalent.backward = pagemtimes(receive.beamformer, channel.backward);
		channel.equivalent.aggregate = channel_aggregate(channel.equivalent.direct, channel.equivalent.forward, channel.equivalent.backward, ris.scatter);
		channel.interference(b, r) = norm(channel.equivalent.aggregate(:, :, ~logical(eye(network.link))), 'fro') ^ 2;
		while ~iter.converge
			iter.interference = channel.interference(b, r);
			receive.beamformer = receive_min_interference(channel.aggregate, transmit.beamformer);
			transmit.beamformer = receive_min_interference(pagectranspose(channel.aggregate), pagectranspose(receive.beamformer));
			[ris.scatter, channel.aggregate] = ris_min_interference();
			channel.interference(b, r) = norm(channel.equivalent.aggregate(:, :, ~logical(eye(network.link))), 'fro') ^ 2;
			iter.converge = (abs(channel.interference(b, r) - iter.interference) <= iter.tolerance);
			iter.counter = iter.counter + 1;
		end
	end
end

channel.power.direct = mean(channel.power.direct, 2);
channel.power.cascaded = mean(channel.power.cascaded, 2);
channel.power.aggregate = mean(channel.power.aggregate, 2);

figure('Name', 'Channel Power vs RIS Group Size', 'Position', [0, 0, 500, 400]);
bar(1, [channel.power.direct, channel.power.cascaded]', 'stacked');
set(gca, 'XTickLabel', []);
hold all;
for b = 1 : number.bond
	bar(b + 1, channel.power.aggregate(b));
end
ylabel('Channel Power [J]');
set(gca, 'XTick', 1 : number.bond + 1);
handle.legend = legend(['Direct', 'Cascaded', '$N_g = ' + string(ris.bond) + '$'], 'Orientation', 'horizontal', 'Location', 'northoutside');
handle.legend.ItemTokenSize = [10, 10];
savefig('plots/max_power_bond.fig');
