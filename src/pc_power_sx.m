clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(8, 2 .^ (2 : 2 : 8), 4);
% [channel.pathloss.distance.direct, channel.pathloss.distance.forward, channel.pathloss.distance.backward, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(-14.7, -10, -6.3, -3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.antenna, number.realization] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 2);

for r = 1 : number.realization
	% * No RIS
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		clear scatter_power_max_pc;
		for b = 1 : number.bond(a)
			reflect.bond = 2 .^ (b - 1);
			[reflect.beamformer, channel.aggregate] = scatter_power_max_pc(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate(b, a, r) = norm(channel.aggregate, 'fro') ^ 2;
		end
	end
end
channel.power.direct = mean(channel.power.direct, ndims(channel.power.direct));
channel.power.aggregate = mean(channel.power.aggregate, ndims(channel.power.aggregate));
save('data/pc_power_sx.mat');

figure('Name', 'Channel Power vs RIS Configuration', 'Position', [0, 0, 500, 400]);
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$');
hold all;
handle.power.direct = refline(0, channel.power.direct);
set(handle.power.direct, 'Color', 'k', 'Marker', 'none', 'DisplayName', '$N^\mathrm{S} = 0$');
for a = 1 : number.antenna
	handle.power.aggregate(a) = plot(1 : number.bond(a), channel.power.aggregate(1 : number.bond(a), a), 'DisplayName', '$N^\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
end
style_plot(handle.power.aggregate);
hold off; grid on; box on; legend('Location', 'nw');
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
savefig('plots/pc_power_sx.fig');
matlab2tikz('../assets/simulation/pc_power_sx.tex', 'width', '8cm', 'height', '6cm');
