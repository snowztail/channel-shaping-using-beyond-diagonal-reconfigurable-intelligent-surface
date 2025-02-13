clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2, 64, 2);
[channel.rank, reflect.bond, channel.uncertainty] = deal(min(transmit.antenna, receive.antenna), [1, reflect.antenna], [1e-2, 1e-1, 5e-1]);
assert(channel.rank == 2, 'The plot function is only for 2-sv case.');
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[channel.weight(1, :, 1), channel.weight(2, :, 1)] = pol2cart(pi / 4 : pi / 64 : pi * 3 / 4, 1);
[channel.weight(1, :, 2), channel.weight(2, :, 2)] = pol2cart(pi / 4 : -pi / 64 : -pi / 4, 1);
[channel.weight(1, :, 3), channel.weight(2, :, 3)] = pol2cart(5 * pi / 4 : pi / 64 : pi * 7 / 4, 1);
[channel.weight(1, :, 4), channel.weight(2, :, 4)] = pol2cart(5 * pi / 4 : -pi / 64 : pi * 3 / 4, 1);
[number.weight, number.case, number.bond, number.uncertainty, number.realization] = deal(size(channel.weight, 2), size(channel.weight, 3), 2, length(channel.uncertainty), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	% * No RIS
	channel.singular.direct(:, r) = svd(channel.direct);
	% * Have RIS
	channel.forward.actual = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward.actual = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	for u = 1 : number.uncertainty
		channel.forward.estimate = estimation_effect(channel.forward.actual, channel.uncertainty(u));
		channel.backward.estimate = estimation_effect(channel.backward.actual, channel.uncertainty(u));
		for b = 1 : number.bond
			for c = 1 : number.case
				clear scatter_singular;
				for w = 1 : number.weight
					[reflect.beamformer] = scatter_singular(channel.direct, channel.forward.estimate, channel.backward.estimate, channel.weight(:, w, c), reflect.bond(b));
					channel.aggregate.actual = channel_aggregate(channel.direct, channel.forward.actual, channel.backward.actual, reflect.beamformer);
					channel.singular.aggregate(:, w, c, b, u, r) = svd(channel.aggregate.actual);
				end
			end
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
channel.singular.aggregate = mean(channel.singular.aggregate, ndims(channel.singular.aggregate));
save('data/pc_singular_csi.mat');

figure('Name', 'Channel Singular Value Pareto Front vs Channel Estimation Error', 'Position', [0, 0, 500, 400]);
handle.color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
hold all;
handle.singular.direct = scatter(channel.singular.direct(1), channel.singular.direct(2), 200, 'k', '^', 'DisplayName', 'Direct');
for u = 1 : number.uncertainty
	point.all = [vec(channel.singular.aggregate(1, :, :, 1, u)), vec(channel.singular.aggregate(2, :, :, 1, u))];
	% * Pareto frontier
	point.pareto = point.all(convhull(point.all), :);
	handle.singular.pareto(1, u) = plot(point.pareto(:, 1), point.pareto(:, 2), 'DisplayName', 'D: $\epsilon = ' + string(channel.uncertainty(u)) + '$', 'LineStyle', '--', 'Color', [handle.color{1} (number.uncertainty - u + 1) / number.uncertainty]);
	% * channel power gain-optimal point
	[~, index.power] = max(vecnorm(point.pareto, 2, 2));
	point.power = point.pareto(index.power, :);
	handle.singular.power(1, u) = scatter(point.power(1), point.power(2), 200, 'MarkerEdgeColor', handle.color{1}, 'MarkerEdgeAlpha', (number.uncertainty - u + 1) / number.uncertainty);
	% * power transfer-optimal point
	[~, index.energy] = max(point.pareto(:, 1));
	point.energy = point.pareto(index.energy, :);
	handle.singular.energy(1, u) = scatter(point.energy(1), point.energy(2), 200, 'Marker', 'x', 'MarkerEdgeColor', handle.color{1}, 'MarkerEdgeAlpha', (number.uncertainty - u + 1) / number.uncertainty);
	% * rate-optimal arc
	[~, index.rate.low] = max(point.pareto(:, 1));
	[~, index.rate.high] = max(point.pareto(:, 2));
	index.rate.range = circshift(setdiff(1 : size(point.pareto, 1), min(index.rate.low, index.rate.high) + 1 : max(index.rate.low, index.rate.high) - 1), -index.rate.high);
	handle.singular.rate(1, u) = plot(point.pareto(index.rate.range, 1), point.pareto(index.rate.range, 2), 'LineStyle', '-', 'Color', [handle.color{1} (number.uncertainty - u + 1) / number.uncertainty]);
end
for u = 1 : number.uncertainty
	point.all = [vec(channel.singular.aggregate(1, :, :, 2, u)), vec(channel.singular.aggregate(2, :, :, 2, u))];
	% * Pareto frontier
	point.pareto = point.all(convhull(point.all), :);
	handle.singular.pareto(2, u) = plot(point.pareto(:, 1), point.pareto(:, 2), 'DisplayName', 'BD: $\epsilon = ' + string(channel.uncertainty(u)) + '$', 'LineStyle', '--', 'Color', [handle.color{2} (number.uncertainty - u + 1) / number.uncertainty]);
	% * channel power gain-optimal point
	[~, index.power] = max(vecnorm(point.pareto, 2, 2));
	point.power = point.pareto(index.power, :);
	handle.singular.power(2, u) = scatter(point.power(1), point.power(2), 200, 'MarkerEdgeColor', handle.color{2}, 'MarkerEdgeAlpha', (number.uncertainty - u + 1) / number.uncertainty);
	% * power transfer-optimal point
	[~, index.energy] = max(point.pareto(:, 1));
	point.energy = point.pareto(index.energy, :);
	handle.singular.energy(2, u) = scatter(point.energy(1), point.energy(2), 200, 'Marker', 'x', 'MarkerEdgeColor', handle.color{2}, 'MarkerEdgeAlpha', (number.uncertainty - u + 1) / number.uncertainty);
	% * rate-optimal arc
	[~, index.rate.low] = max(point.pareto(:, 1));
	[~, index.rate.high] = max(point.pareto(:, 2));
	index.rate.range = circshift(setdiff(1 : size(point.pareto, 1), min(index.rate.low, index.rate.high) + 1 : max(index.rate.low, index.rate.high) - 1), -index.rate.high);
	handle.singular.rate(2, u) = plot(point.pareto(index.rate.range, 1), point.pareto(index.rate.range, 2), 'LineStyle', '-', 'Color', [handle.color{2} (number.uncertainty - u + 1) / number.uncertainty]);
end
handle.dummy.power = scatter(nan, nan, 200, 'MarkerEdgeColor', '#808080', 'DisplayName', 'P-max');
handle.dummy.energy = scatter(nan, nan, 200, 'Marker', 'x', 'MarkerEdgeColor', '#808080', 'DisplayName', 'E-max');
handle.dummy.rate = plot(nan, nan, 'Color', '#808080', 'DisplayName', 'R-max');
hold off; grid on; xlim tight; ylim tight; box on; legend([handle.singular.direct, handle.singular.pareto(1, :), handle.singular.pareto(2, :), handle.dummy.power, handle.dummy.energy, handle.dummy.rate], 'Location', 'nw');
xlabel('$\sigma_1(\mathbf{H})$');
ylabel('$\sigma_2(\mathbf{H})$');
savefig('plots/pc_singular_csi.fig');
matlab2tikz('../assets/simulation/pc_singular_csi.tex', 'width', '8cm', 'height', '6cm', 'extraaxisoptions', {'every axis plot/.append style={line width=1.5pt}'});


function [H] = estimation_effect(H, sigma)
	H = H .* (1 + sigma * sqrt(0.5) * (randn(size(H)) + 1i * randn(size(H))));
end
