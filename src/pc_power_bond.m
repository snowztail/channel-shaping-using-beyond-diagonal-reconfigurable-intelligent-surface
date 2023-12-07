clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(8, 256, 4);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 2 : log2(reflect.antenna)));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.realization, flag.direct] = deal(length(reflect.bond), 1e1, true);

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
for b = 1 : number.bond
	handle.power(b) = bar(b, channel.power.aggregate(b));
end
handle.power(b + 1) = refline(0, channel.power.direct);
handle.power(b + 2) = refline(0, channel.power.cascaded);
set(handle.power(b + 1), 'Color', 'k', 'Marker', 'none');
set(handle.power(b + 2), 'Color', '#ED8198', 'Marker', 'none');
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
legend(['$L = ' + string(vec(reflect.bond)) + '$'; 'Direct'; 'Cascaded'], 'Location', 'nw');
if flag.direct == true
	savefig('plots/pc_power_bond_hd.fig');
else
	savefig('plots/pc_power_bond_nd.fig');
end
