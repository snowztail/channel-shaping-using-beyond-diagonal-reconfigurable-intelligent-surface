clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 32, 4);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), [1, reflect.antenna]);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(1, db2pow(-17.5), db2pow(-17.5));
channel.weight = eye(channel.rank);
[number.weight, number.bond, number.realization] = deal(size(channel.weight, 2), length(reflect.bond), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_los(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.auxiliary = channel_auxiliary(channel.direct, channel.forward);
	channel.singular.direct(:, r) = svd(channel.direct);
	channel.singular.auxiliary(:, r) = svd(channel.auxiliary);
	for w = 1 : number.weight
		for b = 1 : number.bond
			clear scatter_singular_pc;
			[reflect.beamformer.min, channel.aggregate.min] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, -channel.weight(:, w), reflect.bond(b));
			clear scatter_singular_pc;
			[reflect.beamformer.max, channel.aggregate.max] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w), reflect.bond(b));
			channel.singular.aggregate.min(:, w, b, r) = svd(channel.aggregate.min);
			channel.singular.aggregate.max(:, w, b, r) = svd(channel.aggregate.max);
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
channel.singular.auxiliary = mean(channel.singular.auxiliary, ndims(channel.singular.auxiliary));
channel.singular.aggregate.min = mean(channel.singular.aggregate.min, ndims(channel.singular.aggregate.min));
channel.singular.aggregate.max = mean(channel.singular.aggregate.max, ndims(channel.singular.aggregate.max));
save('data/pc_singular_bound.mat');

figure('Name', 'Channel Singular Value vs RIS Group Size', 'Position', [0, 0, 500, 400]);
set(gca, 'XLim', [0, channel.rank + 1], 'XTick', 1 : channel.rank, 'XTickLabel', '$\sigma_' + string(vec(1 : channel.rank)) + '$');
hold all;
handle.singular = bar([channel.singular.direct, pagediag(channel.singular.aggregate.min), pagediag(channel.singular.aggregate.max)], 'FaceColor', 'flat');
handle.singular(1).FaceColor = '#000000';
handle.singular(2).FaceColor = '#94C4F5';
handle.singular(3).FaceColor = '#4260AA';
handle.singular(4).FaceColor = '#F26161';
handle.singular(5).FaceColor = '#EA0909';
set(handle.singular, {'DisplayName'}, cellstr(['Direct'; 'Min, $L = ' + string(vec(reflect.bond)) + '$'; 'Max, $L = ' + string(vec(reflect.bond)) + '$']));

for w = 1 : number.weight
	handle.bound(w) = refline(0, channel.singular.auxiliary(w));
end
set(handle.bound(1), 'Color', '#77AC30', 'Marker', 'none', 'LineStyle', '-');
set(handle.bound(2), 'Color', '#77AC30', 'Marker', 'none', 'LineStyle', '--');
set(handle.bound(3), 'Color', '#77AC30', 'Marker', 'none', 'LineStyle', ':');
set(handle.bound(4), 'Color', '#77AC30', 'Marker', 'none', 'LineStyle', '-.');
set(handle.bound, {'DisplayName'}, cellstr('$\sigma_' + string(vec(1 : channel.rank)) + '(\mathbf{T})$'));
hold off; grid on; box on; legend('Location', 'ne', 'NumColumns', 2);
ylabel('Amplitude');
savefig('plots/pc_singular_los.fig');
matlab2tikz('../assets/simulation/pc_singular_bound.tex', 'width', '10cm', 'height', '7.5cm', 'extraaxisoptions', {'legend columns=5', 'transpose legend', 'legend style={/tikz/column 2/.style={column sep=5pt}}'});

function [H_a] = channel_auxiliary(H_d, H_f)
	[~, ~, v_f] = svds(H_f, 1);
	H_a = sqrtm(H_d * (eye(size(H_d, 2)) - v_f * v_f') * H_d');
end

function [X] = pagediag(A)
	[M, ~, N] = size(A);
	X = zeros(M, N);
	for n = 1 : N
		X(:, n) = diag(A(:, :, n));
	end
end
