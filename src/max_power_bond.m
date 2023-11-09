clc; clear; close; setup;

[base.antenna, ris.antenna, user.antenna] = deal(4, 256, 2);
ris.bond = 2 .^ (0 : 2 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.rank = min(base.antenna, user.antenna);
channel.weight = simplex_standard(channel.rank, 0.1);
[number.bond, number.weight, number.realization] = deal(length(ris.bond), length(channel.weight), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_ricean(base.antenna, 'ula', user.antenna, 'ula');
	% channel.direct = 0 * channel.direct;
	channel.forward = sqrt(channel.pathloss.forward) * fading_ricean(base.antenna, 'ula', ris.antenna, 'upa');
	channel.backward = sqrt(channel.pathloss.backward) * fading_ricean(ris.antenna, 'upa', user.antenna, 'ula');
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	channel.power.cascaded(r) = sum(svds(channel.forward, channel.rank) .^ 2 .* svds(channel.backward, channel.rank) .^ 2);
	for b = 1 : number.bond
		ris.scatter = eye(ris.antenna);
		channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
		[ris.scatter, channel.aggregate] = ris_max_power(channel.direct, channel.forward, channel.backward, ris.scatter, ris.group(b));
		channel.power.aggregate(b, r) = norm(channel.aggregate, 'fro') ^ 2;
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
