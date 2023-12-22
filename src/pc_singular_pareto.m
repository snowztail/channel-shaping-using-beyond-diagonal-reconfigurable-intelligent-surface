clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 64, 2);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 2 .^ (0 : 2 : log2(reflect.antenna)));
assert(channel.rank == 2, 'The plot function is only for 2-sv case.');
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(1, db2pow(-17.5), db2pow(-17.5));
channel.weight = simplex_standard(channel.rank, 0.1);
[number.weight, number.bond, number.realization] = deal(size(channel.weight, 2), length(reflect.bond), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.singular.direct(:, r) = svd(channel.direct);
	for b = 1 : number.bond
		clear scatter_singular_min_pc scatter_singular_max_pc;
		for w = 1 : number.weight
			[reflect.beamformer.min, channel.aggregate.min] = scatter_singular_min_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w), reflect.bond(b));
			[reflect.beamformer.max, channel.aggregate.max] = scatter_singular_max_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w), reflect.bond(b));
			channel.singular.aggregate.min(:, w, b, r) = svd(channel.aggregate.min);
			channel.singular.aggregate.max(:, w, b, r) = svd(channel.aggregate.max);
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
channel.singular.aggregate.min = mean(channel.singular.aggregate.min, ndims(channel.singular.aggregate.min));
channel.singular.aggregate.max = mean(channel.singular.aggregate.max, ndims(channel.singular.aggregate.max));
save('data/pc_singular_pareto.mat');

figure('Name', 'Channel Singular Value Pareto Front vs RIS Group Size', 'Position', [0, 0, 500, 400]);
hold all;
handle.singular.direct = scatter(channel.singular.direct(1), channel.singular.direct(2), 200, 'k', '^', 'DisplayName', 'Direct');
for b = 1 : number.bond
	handle.singular.pareto.min(b) = plot(channel.singular.aggregate.min(1, :, b), channel.singular.aggregate.min(2, :, b), 'DisplayName', '$L = ' + string(reflect.bond(b)) + '$');
	handle.singular.pareto.max(b) = plot(channel.singular.aggregate.max(1, :, b), channel.singular.aggregate.max(2, :, b), 'DisplayName', '$L = ' + string(reflect.bond(b)) + '$');
end
style_plot(handle.singular.pareto.min);
style_plot(handle.singular.pareto.max);
hold off; grid on; xlim tight; ylim tight; box on; legend([handle.singular.direct, handle.singular.pareto.min], 'Location', 'nw');
xlabel('$\sigma_1$');
ylabel('$\sigma_2$');
savefig('plots/pc_singular_pareto.fig');
matlab2tikz('../assets/simulation/pc_singular_pareto.tex', 'width', '8cm', 'height', '6cm', 'extraaxisoptions', {'every axis plot/.append style={line width=1.5pt}'});
