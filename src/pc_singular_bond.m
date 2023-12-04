clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna] = deal(4, 256, 2);
ris.bond = 2 .^ (0 : 2 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.rank = min(transmit.antenna, receive.antenna);
channel.weight = simplex_standard(channel.rank, 0.1);
[number.bond, number.weight, number.realization] = deal(length(ris.bond), length(channel.weight), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_ricean(receive.antenna, 'ula', transmit.antenna, 'ula');
	% channel.direct = 0 * channel.direct;
	channel.forward = sqrt(channel.pathloss.forward) * fading_ricean(ris.antenna, 'upa', transmit.antenna, 'ula');
	channel.backward = sqrt(channel.pathloss.backward) * fading_ricean(receive.antenna, 'ula', ris.antenna, 'upa');
	channel.sv.direct(:, r) = svd(channel.direct);
	for w = 1 : number.weight
		for b = 1 : number.bond
			ris.scatter = eye(ris.antenna);
			channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
			[ris.scatter, channel.aggregate] = ris_max_wsv(channel.direct, channel.forward, channel.backward, channel.weight(:, w), ris.scatter, ris.group(b));
			channel.sv.aggregate(:, b, w, r) = svd(channel.aggregate);
		end
	end
end

channel.sv.direct = mean(channel.sv.direct, 2);
channel.sv.aggregate = mean(channel.sv.aggregate, 4);

figure('Name', 'Channel Singular Value vs RIS Group Size', 'Position', [0, 0, 500, 400]);
number.plot = 3;
handle.window = tiledlayout(number.plot, 1);
for p = 1 : number.plot
	w = (p - 1) * (number.weight - 1) / (number.plot - 1) + 1;
	handle.axis(p) = nexttile;
	bar([channel.sv.direct, channel.sv.aggregate(:, :, w)]);
	title('$\rho_1 = ' + string(channel.weight(1, w)) + '$, $\rho_2 = ' + string(channel.weight(2, w)) + '$');
end
handle.xlabel = xlabel(handle.window, 'Index');
handle.ylabel = ylabel(handle.window, 'Singular Value');
handle.xlabel.Interpreter= 'latex'; handle.ylabel.Interpreter= 'latex';
handle.legend = legend(['No RIS', '$N_g = ' + string(ris.bond) + '$'], 'Orientation', 'horizontal');
handle.legend.ItemTokenSize = [10, 10]; handle.legend.Layout.Tile = 'north';
linkaxes(handle.axis);
savefig('plots/max_wsv_bond.fig');
