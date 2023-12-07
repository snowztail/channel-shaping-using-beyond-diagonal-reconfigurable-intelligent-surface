clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(8, 2 .^ [-inf, 2 : 2 : 8], 4);
% [channel.pathloss.distance.direct, channel.pathloss.distance.forward, channel.pathloss.distance.backward, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(-14.7, -10, -6.3, -3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.antenna, number.realization] = deal(length(reflect.antenna), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for a = 1 : number.antenna
		if reflect.antenna(a) == 0
			channel.power.aggregate(1, a, r) = norm(channel.direct, 'fro') ^ 2;
		else
			channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
			channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
			reflect.bond = 2 .^ (0 : log2(reflect.antenna(a))); number.bond = length(reflect.bond);
			for b = 1 : number.bond
				[reflect.beamformer, channel.aggregate] = scatter_power_pc(channel.direct, channel.forward, channel.backward, eye(reflect.antenna(a)), reflect.bond(b));
				channel.power.aggregate(b, a, r) = norm(channel.aggregate, 'fro') ^ 2;
			end
		end
	end
end
channel.power.aggregate = mean(channel.power.aggregate, 3);
save('data/pc_power_sx.mat');

figure('Name', 'Channel Power vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
set(gca, 'XLim', [0, log2(max(reflect.antenna))], 'XTick', 0 : log2(max(reflect.antenna)), 'XTickLabel', '$2^' + string(vec(0 : log2(max(reflect.antenna)))) + '$');
for a = 1 : number.antenna
	if reflect.antenna(a) == 0
		handle.power(a) = refline(0, channel.power.aggregate(1, a));
	else
		reflect.bond = 2 .^ (0 : log2(reflect.antenna(a))); number.bond = length(reflect.bond);
		handle.power(a) = plot(log2(reflect.bond), channel.power.aggregate(1 : number.bond, a));
	end
end
hold off; grid on;
style_plot(handle.power); set(handle.power(find(reflect.antenna == 0)), 'Color', 'k', 'Marker', 'none');
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
legend('$N_s = ' + string(vec(reflect.antenna)) + '$', 'Location', 'nw');
savefig('plots/pc_power_sx.fig');
