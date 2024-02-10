clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2, 32, 2);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 1 : log2(reflect.antenna)));
assert(channel.rank == 2, 'The plot function is only for 2-sv case.');
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[channel.weight(1, :, 1), channel.weight(2, :, 1)] = pol2cart(pi / 4 : pi / 64 : pi * 3 / 4, 1);
[channel.weight(1, :, 2), channel.weight(2, :, 2)] = pol2cart(pi / 4 : -pi / 64 : -pi / 4, 1);
[channel.weight(1, :, 3), channel.weight(2, :, 3)] = pol2cart(5 * pi / 4 : pi / 64 : pi * 7 / 4, 1);
[channel.weight(1, :, 4), channel.weight(2, :, 4)] = pol2cart(5 * pi / 4 : -pi / 64 : pi * 3 / 4, 1);
[number.weight, number.case, number.bond, number.realization] = deal(size(channel.weight, 2), size(channel.weight, 3), length(reflect.bond), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.singular.direct(:, r) = svd(channel.direct);
	for b = 1 : number.bond
		for c = 1 : number.case
			clear scatter_singular_pc;
			for w = 1 : number.weight
				[reflect.beamformer, channel.aggregate] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w, c), reflect.bond(b));
				channel.singular.aggregate(:, w, c, b, r) = svd(channel.aggregate);
			end
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
channel.singular.aggregate = mean(channel.singular.aggregate, ndims(channel.singular.aggregate));
save('data/pc_singular_pareto.mat');

figure('Name', 'Channel Singular Value Pareto Front vs RIS Group Size', 'Position', [0, 0, 500, 400]);
hold all;
handle.singular.direct = scatter(channel.singular.direct(1), channel.singular.direct(2), 200, 'k', '^', 'DisplayName', 'Direct');
for b = 1 : number.bond
	p = [vec(channel.singular.aggregate(1, :, :, b)), vec(channel.singular.aggregate(2, :, :, b))];
	k = convhull(p);
	handle.singular.pareto(b) = plot(p(k, 1), p(k, 2), 'DisplayName', '$L = ' + string(reflect.bond(b)) + '$');
end
style_plot(handle.singular.pareto);
hold off; grid on; xlim tight; ylim tight; box on; legend([handle.singular.direct, handle.singular.pareto], 'Location', 'nw');
xlabel('$\sigma_1$');
ylabel('$\sigma_2$');
savefig('plots/pc_singular_pareto.fig');
matlab2tikz('../assets/simulation/pc_singular_pareto.tex', 'width', '8cm', 'height', '6cm', 'extraaxisoptions', {'every axis plot/.append style={line width=1.5pt}'});
