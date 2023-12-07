clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna] = deal(4, 16, 2);
[channel.rank, ris.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 2 : log2(ris.antenna)));
assert(channel.rank == 2, 'The plot function is only for 2-sv case.');
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(1, 0.1, 0.1);
channel.weight = simplex_standard(channel.rank, 0.1);
[number.weight, number.bond, number.realization] = deal(size(channel.weight, 2), length(ris.bond), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(ris.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, ris.antenna);
	channel.singular.direct(:, r) = svd(channel.direct);
	for b = 1 : number.bond
		ris.scatter = scatter_power_pc(channel.direct, channel.forward, channel.backward, eye(ris.antenna), ris.bond(b));
		for w = 1 : number.weight
			channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
			[ris.scatter, channel.aggregate] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w), ris.scatter, ris.bond(b));
			channel.singular.aggregate(:, w, b, r) = svd(channel.aggregate);
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, 2);
channel.singular.aggregate = mean(channel.singular.aggregate, 4);
save('data/pc_singular_bond.mat');

figure('Name', 'Channel Singular Value vs RIS Group Size', 'Position', [0, 0, 500, 400]);
hold all;
for s = 1 : channel.rank
	handle.singular(1, s) = refline(0, channel.singular.direct(s));
	handle.singular(1, s).DisplayName = '$\sigma_' + string(s) + '(\mathbf{H}^\mathrm{D})$';
	for b = 1 : number.bond
		handle.singular(b + 1, s) = plot(channel.weight(1, :), channel.singular.aggregate(s, :, b), 'DisplayName', '$\sigma_' + string(s) + '(\mathbf{H})$, $L = ' + string(ris.bond(b)) + '$');
	end
end
hold off; grid on; ylim auto;
style_plot(handle.singular, number.bond + 1);
xlabel('$\rho_1$');
ylabel('Singular Value');
legend('Location', 'sw');
savefig('plots/pc_singular_bond.fig');
