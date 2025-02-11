clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2, 32, 2);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), unique([1, 4, 16, reflect.antenna]));
[transmit.power, receive.noise] = deal(db2pow([-20, 20]), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
assert(channel.rank == 2, 'The plot function is only for 2-sv case.');
[channel.weight(1, :, 1), channel.weight(2, :, 1)] = pol2cart(pi / 4 : pi / 256 : pi * 3 / 4, 1);
[channel.weight(1, :, 2), channel.weight(2, :, 2)] = pol2cart(pi / 4 : -pi / 256 : -pi / 4, 1);
[channel.weight(1, :, 3), channel.weight(2, :, 3)] = pol2cart(5 * pi / 4 : pi / 256 : pi * 7 / 4, 1);
[channel.weight(1, :, 4), channel.weight(2, :, 4)] = pol2cart(5 * pi / 4 : -pi / 256 : pi * 3 / 4, 1);
[number.weight, number.case, number.bond, number.power, number.realization, flag.direct] = deal(size(channel.weight, 2), size(channel.weight, 3), length(reflect.bond), length(transmit.power), 1, false);

for r = 1 : number.realization
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna) + eps;
	channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.singular.direct(:, r) = svd(channel.direct);
	for b = 1 : number.bond
		% * Pareto frontier
		for c = 1 : number.case
			for w = 1 : number.weight
				clear scatter_singular;
				[~, channel.aggregate] = scatter_singular(channel.direct, channel.forward, channel.backward, channel.weight(:, w, c), reflect.bond(b));
				channel.singular.aggregate(:, w, c, b, r) = svd(channel.aggregate);
				for p = 1 : number.power
					transmit.beamformer = precoder_rate(channel.aggregate, transmit.power(p), receive.noise);
					receive.rate(p, w, c, b, r) = rate_mimo(channel.aggregate, transmit.beamformer, receive.noise);
				end
			end
		end
		% * rate-optimal arc
		for p = 1 : number.power
			[~, i] = max(vec(receive.rate(p, :, :, b, r)));
			[w, c] = ind2sub([number.weight, number.case], i);
			channel.singular.rate(:, p, b, r) = channel.singular.aggregate(:, w, c, b, r);
		end
		% * power-optimal point
		clear scatter_power_max;
		[~, channel.aggregate] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond(b));
		channel.singular.power(:, b, r) = svd(channel.aggregate);
	end
end
% channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
% channel.singular.aggregate = mean(channel.singular.aggregate, ndims(channel.singular.aggregate));
% channel.singular.power = mean(channel.singular.power, ndims(channel.singular.power));
% channel.singular.rate = mean(channel.singular.rate, ndims(channel.singular.rate));
save('data/pc_singular_pareto.mat');

figure('Name', 'Channel Singular Value Pareto Front vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.color = {'#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEE'; '#A2142F'};
hold all;
handle.singular.direct = scatter(channel.singular.direct(1), channel.singular.direct(2), 200, 'k', '^', 'DisplayName', 'Direct');
for b = 1 : number.bond
	p = [vec(channel.singular.aggregate(1, :, :, b)), vec(channel.singular.aggregate(2, :, :, b))];
	k = convhull(p);
	handle.singular.pareto(b) = plot(p(k, 1), p(k, 2), 'DisplayName', '$L = ' + string(reflect.bond(b)) + '$', 'LineStyle', ':', 'Color', handle.color{b});
	handle.singular.power(b) = scatter(channel.singular.power(1, b), channel.singular.power(2, b), 'MarkerEdgeColor', handle.color{b});
	handle.singular.rate(b) = plot(squeeze(channel.singular.rate(1, :, b)), squeeze(channel.singular.rate(2, :, b)), 'LineStyle', '-', 'Color', handle.color{b});
end
% style_plot(handle.singular.pareto);
hold off; grid on; xlim tight; ylim tight; box on; legend([handle.singular.direct, handle.singular.pareto], 'Location', 'nw');
xlabel('$\sigma_1(\mathbf{H})$');
ylabel('$\sigma_2(\mathbf{H})$');
savefig('plots/pc_singular_pareto.fig');
matlab2tikz('../assets/simulation/pc_singular_pareto.tex', 'width', '8cm', 'height', '6cm', 'extraaxisoptions', {'every axis plot/.append style={line width=1.5pt}'});
