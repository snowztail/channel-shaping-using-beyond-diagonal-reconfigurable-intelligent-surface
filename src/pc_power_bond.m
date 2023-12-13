clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(8, 256, 4);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 2 : log2(reflect.antenna)));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.realization, flag.direct] = deal(length(reflect.bond), 1e1, false);

for r = 1 : number.realization
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	channel.power.cascaded(r) = sum(svds(channel.backward, channel.rank) .^ 2 .* svds(channel.forward, channel.rank) .^ 2);
	for b = 1 : number.bond
		[reflect.beamformer, channel.aggregate] = scatter_power_pc(channel.direct, channel.forward, channel.backward, eye(reflect.antenna), reflect.bond(b));
		channel.power.aggregate(b, r) = norm(channel.aggregate, 'fro') ^ 2;
	end
end
channel.power.direct = mean(channel.power.direct, 2);
channel.power.cascaded = mean(channel.power.cascaded, 2);
channel.power.aggregate = mean(channel.power.aggregate, 2);
if flag.direct == true
	save('data/pc_power_bond_hd.mat');
else
	save('data/pc_power_bond_nd.mat');
end

figure('Name', 'Channel Power vs RIS Group Size', 'Position', [0, 0, 500, 400]);
hold all;
handle.power(1) = bar(channel.power.aggregate, 'FaceColor', '#77AC30');
handle.power(2) = refline(0, channel.power.direct);
handle.power(3) = refline(0, channel.power.cascaded);
hold off; grid on; box on;
set(handle.power(2), 'Color', 'k', 'Marker', 'none');
set(handle.power(3), 'Color', '#ED8198', 'Marker', 'none');
xlim([0, number.bond + 1]);
xticks(1 : number.bond);
xticklabels(reflect.bond);
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
legend(handle.power, {'Composite'; 'Direct'; 'Cascaded'}, 'Location', 'nw');
if flag.direct == true
	savefig('plots/pc_power_bond_hd.fig');
	matlab2tikz('../assets/simulation/pc_power_bond_hd.tex', 'width', '8cm', 'height', '6cm');
else
	savefig('plots/pc_power_bond_nd.fig');
	matlab2tikz('../assets/simulation/pc_power_bond_nd.tex', 'width', '8cm', 'height', '6cm');
end
