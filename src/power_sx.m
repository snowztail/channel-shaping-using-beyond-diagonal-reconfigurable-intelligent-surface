clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(16, 2 .^ (4 : 8), 16);
% [channel.pathloss.distance.direct, channel.pathloss.distance.forward, channel.pathloss.distance.backward, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(-14.7, -10, -6.3, -3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.antenna, number.realization, flag.direct] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 1e1, true);

for r = 1 : number.realization
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		clear scatter_power_max;
		for b = 1 : number.bond(a)
			reflect.bond = 2 .^ (b - 1);
			[reflect.beamformer, channel.aggregate] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate(b, a, r) = norm(channel.aggregate, 'fro') ^ 2;
		end
	end
end
channel.power.direct = mean(channel.power.direct, ndims(channel.power.direct));
channel.power.aggregate = mean(channel.power.aggregate, ndims(channel.power.aggregate));
save('data/power_sx.mat');

figure('Name', 'Channel Power vs RIS Configuration', 'Position', [0, 0, 500, 400]);
if flag.direct
    handle.power.direct = semilogy(1 : max(number.bond), repmat(channel.power.direct, [1, max(number.bond)]), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
    hold on;
end
for a = 1 : number.antenna
	handle.power.aggregate(a) = semilogy(1 : number.bond(a), channel.power.aggregate(1 : number.bond(a), a), 'DisplayName', '$N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
end
style_plot(handle.power.aggregate);
hold off; grid on; box on; legend('Location', 'nw');
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$', 'YLim', [channel.power.direct * 0.95, max(vec(channel.power.aggregate)) * 1.05]);
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
savefig('plots/power_sx.fig');
matlab2tikz('../assets/simulation/power_sx.tex', 'width', '8cm', 'height', '6cm');
