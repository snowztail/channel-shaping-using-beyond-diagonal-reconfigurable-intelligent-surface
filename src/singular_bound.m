clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 256, 4);
[channel.rank.direct, channel.rank.forward, channel.rank.backward, reflect.bond] = deal(min(transmit.antenna, receive.antenna), 4, min(reflect.antenna, receive.antenna), [1, reflect.antenna]);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.weight = eye(channel.rank.direct);
[number.weight, number.bond, number.realization] = deal(size(channel.weight, 2), length(reflect.bond), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_rankn(reflect.antenna, transmit.antenna, channel.rank.forward);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.auxiliary = channel_auxiliary(channel.direct, channel.forward);
	channel.singular.direct(:, r) = svd(channel.direct);
	channel.singular.auxiliary(:, r) = svd(channel.auxiliary);
	for b = 1 : number.bond
		for w = 1 : number.weight
			clear scatter_singular;
			[reflect.beamformer.min, channel.aggregate.min] = scatter_singular(channel.direct, channel.forward, channel.backward, -channel.weight(:, w), reflect.bond(b));
			clear scatter_singular;
			[reflect.beamformer.max, channel.aggregate.max] = scatter_singular(channel.direct, channel.forward, channel.backward, channel.weight(:, w), reflect.bond(b));
            [s.min, s.max] = deal(svd(channel.aggregate.min), svd(channel.aggregate.max));
			channel.singular.aggregate.min(w, b, r) = s.min(w);
			channel.singular.aggregate.max(w, b, r) = s.max(w);
		end
	end
end
channel.singular.direct = mean(channel.singular.direct, ndims(channel.singular.direct));
channel.singular.auxiliary = mean(channel.singular.auxiliary, ndims(channel.singular.auxiliary));
channel.singular.aggregate.min = mean(channel.singular.aggregate.min, ndims(channel.singular.aggregate.min));
channel.singular.aggregate.max = mean(channel.singular.aggregate.max, ndims(channel.singular.aggregate.max));
save('data/singular_bound.mat');

handle.singular = plotBarStackGroups(cat(3, channel.singular.aggregate.min, channel.singular.direct - channel.singular.aggregate.min, channel.singular.aggregate.max - channel.singular.direct), cellstr('$\sigma_' + string(vec(1 : channel.rank.direct)) + '(\mathbf{H})$'));
for w = 1 : number.weight
	handle.bound(w) = refline(0, channel.singular.auxiliary(w));
end
set(handle.bound, {'Color'}, {'#77AC30'}, {'LineStyle'}, {'-'; '--'; ':'; '-.'}, {'DisplayName'}, cellstr('$\sigma_' + string(vec(1 : channel.rank.direct)) + '(\mathbf{T})$'));
hold off; grid on; box on; legend('Location', 'ne', 'NumColumns', 2);
ylabel('Amplitude');
savefig('plots/singular_bound.fig');
matlab2tikz('../assets/simulation/singular_bound.tex', 'width', '10cm', 'height', '7.5cm', 'extraaxisoptions', {'legend columns=4', 'transpose legend', 'legend style={/tikz/column 2/.style={column sep=5pt}}'});

function [H_a] = channel_auxiliary(H_d, H_f)
	[~, ~, V_f] = svds(H_f, rank(H_f));
	H_a = sqrtm(H_d * (eye(size(H_d, 2)) - V_f * V_f') * H_d');
end

function [h] = plotBarStackGroups(stackData, groupLabels)
	NumGroupsPerAxis = size(stackData, 1);
	NumStacksPerGroup = size(stackData, 2);
	groupBins = 1 : NumGroupsPerAxis;
	MaxGroupWidth = 0.65;
	groupOffset = MaxGroupWidth / NumStacksPerGroup;
	figure
	hold all;
	for i = 1 : NumStacksPerGroup
		Y = squeeze(stackData(:, i, :));
		internalPosCount = i - ((NumStacksPerGroup + 1) / 2);
		groupDrawPos = internalPosCount * groupOffset + groupBins;
		h(i, :) = bar(Y, 'stacked');
		set(h(i, :),'BarWidth',groupOffset);
		set(h(i, :),'XData',groupDrawPos);
	end
	hold off;
	set(gca,'XTickMode','manual');
	set(gca,'XTick',1:NumGroupsPerAxis);
	set(gca,'XTickLabelMode','manual');
	set(gca,'XTickLabel',groupLabels);
	set(h(:, 1),'FaceAlpha',0,'EdgeAlpha',0);
	set(get(get(h(1, 1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	set(get(get(h(2, 1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	set(h(1, 2),'FaceColor','#0072BD','FaceAlpha',0.5,'DisplayName','D-min');
	set(h(2, 2),'FaceColor','#0072BD','DisplayName','BD-min');
	set(h(1, 3),'FaceColor','#D95319','FaceAlpha',0.5,'DisplayName','D-max');
	set(h(2, 3),'FaceColor','#D95319','DisplayName','BD-max');
end
