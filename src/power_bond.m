clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(1, 256, 1);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 2 : log2(reflect.antenna)));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.realization, flag.direct] = deal(length(reflect.bond), 1e2, false);

for r = 1 : number.realization
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	% * Have RIS
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.power.cascaded(r) = sum(svds(channel.backward, channel.rank) .^ 2 .* svds(channel.forward, channel.rank) .^ 2);
	clear scatter_power_max;
	for b = 1 : number.bond
		[reflect.beamformer, channel.aggregate] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond(b));
		channel.power.aggregate(b, r) = norm(channel.aggregate, 'fro') ^ 2;
	end
end
channel.power.direct = mean(channel.power.direct, ndims(channel.power.direct));
channel.power.cascaded = mean(channel.power.cascaded, ndims(channel.power.cascaded));
channel.power.aggregate = mean(channel.power.aggregate, ndims(channel.power.aggregate));
save('data/power_bond.mat');

figure('Name', 'Channel Power vs RIS Group Size', 'Position', [0, 0, 500, 400]);
set(gca, 'XLim', [0, number.bond + 1], 'XTick', 1 : number.bond, 'XTickLabel', reflect.bond, 'YLim', [0, max(channel.power.aggregate) * 1.1]);
hold all;
handle.power.aggregate = bar(channel.power.aggregate, 'FaceColor', '#77AC30', 'DisplayName', 'Equivalent');
if flag.direct
	handle.power.direct = refline(0, channel.power.direct);
	set(handle.power.direct, 'Color', 'k', 'Marker', 'none', 'LineStyle', '-', 'DisplayName', 'Direct');
end
handle.power.cascaded = refline(0, channel.power.cascaded);
set(handle.power.cascaded, 'Color', '#ED8198', 'Marker', 'none', 'LineStyle', '--', 'DisplayName', 'Cascaded');
hold off; grid on; box on; legend('Location', 'se');
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
savefig('plots/power_bond.fig');
matlab2tikz('../assets/simulation/power_bond.tex', 'width', '8cm', 'height', '6cm');
