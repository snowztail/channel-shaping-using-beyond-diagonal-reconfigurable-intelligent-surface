clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(16, 2 .^ (4 : 8), 16);
% [channel.pathloss.distance.direct, channel.pathloss.distance.forward, channel.pathloss.distance.backward, channel.pathloss.exponent.direct, channel.pathloss.exponent.forward, channel.pathloss.exponent.backward] = deal(-14.7, -10, -6.3, -3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.antenna, number.realization, flag.direct] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 1e2, false);

for r = 1 : number.realization
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		for b = 1 : number.bond(a)
			clear scatter_power_max;
			reflect.bond = 2 .^ (b - 1);
			[reflect.beamformer.iterative, channel.aggregate.iterative] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate.iterative(b, a, r) = norm(channel.aggregate.iterative, 'fro') ^ 2;
		end
		if flag.direct
			[reflect.beamformer.procrustes.left, channel.aggregate.procrustes.left] = scatter_procrustes(channel.direct, channel.forward, channel.backward, "left");
			[reflect.beamformer.procrustes.right, channel.aggregate.procrustes.right] = scatter_procrustes(channel.direct, channel.forward, channel.backward, "right");
			channel.power.aggregate.procrustes.left(a, r) = norm(channel.aggregate.procrustes.left, 'fro') ^ 2;
			channel.power.aggregate.procrustes.right(a, r) = norm(channel.aggregate.procrustes.right, 'fro') ^ 2;
		else
			[reflect.beamformer.explicit, channel.aggregate.explicit] = scatter_power_max_explicit(channel.forward, channel.backward);
			channel.power.aggregate.explicit(a, r) = norm(channel.aggregate.explicit, 'fro') ^ 2;
		end
	end
end
channel.power.direct = mean(channel.power.direct, ndims(channel.power.direct));
channel.power.aggregate.iterative = mean(channel.power.aggregate.iterative, ndims(channel.power.aggregate.iterative));
if flag.direct
    channel.power.aggregate.procrustes.left = mean(channel.power.aggregate.procrustes.left, ndims(channel.power.aggregate.procrustes.left));
    channel.power.aggregate.procrustes.right = mean(channel.power.aggregate.procrustes.right, ndims(channel.power.aggregate.procrustes.right));
else
	channel.power.aggregate.explicit = mean(channel.power.aggregate.explicit, ndims(channel.power.aggregate.explicit));
end
% save('data/pc_power_sx.mat');

figure('Name', 'Channel Power vs RIS Configuration', 'Position', [0, 0, 500, 400]);
if flag.direct
    handle.power.direct = semilogy(1 : max(number.bond), repmat(channel.power.direct, [1, max(number.bond)]), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
    hold on;
	for a = 1 : number.antenna
		handle.power.aggregate.procrustes.left(a) = scatter(number.bond(a), channel.power.aggregate.procrustes.left(a), 'Marker', '<');
		hold on;
		handle.power.aggregate.procrustes.right(a) = scatter(number.bond(a), channel.power.aggregate.procrustes.right(a), 'Marker', '>');
		hold on;
	end
	handle.power.aggregate.procrustes.dummy.left = scatter(nan, nan, 'MarkerEdgeColor', '#808080', 'Marker', '<', 'DisplayName', 'OP-left');
	hold on;
	handle.power.aggregate.procrustes.dummy.right = scatter(nan, nan, 'MarkerEdgeColor', '#808080', 'Marker', '>', 'DisplayName', 'OP-right');
	hold on;
	set(handle.power.aggregate.procrustes.left, {'MarkerEdgeColor'}, {'#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'});
	set(handle.power.aggregate.procrustes.right, {'MarkerEdgeColor'}, {'#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'});
else
	for a = 1 : number.antenna
		handle.power.aggregate.explicit(a) = scatter(number.bond(a), channel.power.aggregate.explicit(a), 'Marker', 'o');
		hold on;
	end
	handle.power.aggregate.dummy.explicit = scatter(nan, nan, 'MarkerEdgeColor', '#808080', 'Marker', 'o', 'DisplayName', 'Explicit');
	hold on;
	set(handle.power.aggregate.explicit, {'MarkerEdgeColor'}, {'#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'});
end
for a = 1 : number.antenna
	handle.power.aggregate.iterative(a) = semilogy(1 : number.bond(a), channel.power.aggregate.iterative(1 : number.bond(a), a), 'DisplayName', '$N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
end
style_plot(handle.power.aggregate.iterative);
set(handle.power.aggregate.iterative, {'Marker'}, {'none'});
if flag.direct
	legend([handle.power.direct, handle.power.aggregate.procrustes.dummy.left, handle.power.aggregate.procrustes.dummy.right, handle.power.aggregate.iterative], 'Location', 'nw');
else
	legend([handle.power.aggregate.dummy.explicit, handle.power.aggregate.iterative], 'Location', 'nw');
end
hold off; grid on; box on;
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$', 'YLim', [channel.power.direct * 0.95, max(vec(channel.power.aggregate.iterative)) * 1.05]);
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
% savefig('plots/pc_power_sx.fig');
% matlab2tikz('../assets/simulation/pc_power_sx.tex', 'width', '8cm', 'height', '6cm');


function [Theta, H] = scatter_power_max_explicit(H_f, H_b)
	[U_f, ~, ~] = svd(H_f);
	[~, ~, V_b] = svd(H_b);
	Theta = V_b * U_f';
	H = H_b * Theta * H_f;
end
