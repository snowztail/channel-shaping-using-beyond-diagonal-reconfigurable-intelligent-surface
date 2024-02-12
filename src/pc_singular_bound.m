clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 64, 4);
[channel.rank, reflect.bond] = deal(min(transmit.antenna, receive.antenna), [1, reflect.antenna]);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.weight = eye(channel.rank);
[number.weight, number.bond, number.realization] = deal(size(channel.weight, 2), length(reflect.bond), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.forward = sqrt(channel.pathloss.forward) * fading_los(reflect.antenna, transmit.antenna);
	channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	channel.auxiliary = channel_auxiliary(channel.direct, channel.forward);
	channel.singular.direct(:, r) = svd(channel.direct);
	channel.singular.auxiliary(:, r) = svd(channel.auxiliary);
	for b = 1 : number.bond
		for w = 1 : number.weight
			clear scatter_singular_pc;
			[reflect.beamformer.min, channel.aggregate.min] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, -channel.weight(:, w), reflect.bond(b));
			clear scatter_singular_pc;
			[reflect.beamformer.max, channel.aggregate.max] = scatter_singular_pc(channel.direct, channel.forward, channel.backward, channel.weight(:, w), reflect.bond(b));
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
save('data/pc_singular_bound.mat');

handle.singular = plotBarStackGroups(cat(3, channel.singular.aggregate.min, channel.singular.direct - channel.singular.aggregate.min, channel.singular.aggregate.max - channel.singular.direct), cellstr('$\sigma_' + string(vec(1 : channel.rank)) + '(\mathbf{H})$'));
for w = 1 : number.weight
	handle.bound(w) = refline(0, channel.singular.auxiliary(w));
end
set(handle.bound, {'Color'}, {'#77AC30'}, {'LineStyle'}, {'-'; '--'; ':'; '-.'}, {'DisplayName'}, cellstr('$\sigma_' + string(vec(1 : channel.rank)) + '(\mathbf{T})$'));
hold off; grid on; box on; legend('Location', 'ne', 'NumColumns', 2);
ylabel('Amplitude');
savefig('plots/pc_singular_bound.fig');
matlab2tikz('../assets/simulation/pc_singular_bound.tex', 'width', '10cm', 'height', '7.5cm', 'extraaxisoptions', {'legend columns=4', 'transpose legend', 'legend style={/tikz/column 2/.style={column sep=5pt}}'});

function [H_a] = channel_auxiliary(H_d, H_f)
	[~, ~, v_f] = svds(H_f, 1);
	H_a = sqrtm(H_d * (eye(size(H_d, 2)) - v_f * v_f') * H_d');
end

function [h] = plotBarStackGroups(stackData, groupLabels)
%% Plot a set of stacked bars, but group them according to labels provided.
%%
%% Params: 
%%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%%
%% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
%%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
figure
hold all; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked');
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
end
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
set(h(:, 1),'FaceAlpha',0,'EdgeAlpha',0);
set(get(get(h(1, 1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h(2, 1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(h(1, 2),'FaceColor','#0072BD','FaceAlpha',0.5,'DisplayName','D-max');
set(h(2, 2),'FaceColor','#0072BD','DisplayName','BD-max');
set(h(1, 3),'FaceColor','#D95319','FaceAlpha',0.5,'DisplayName','D-min');
set(h(2, 3),'FaceColor','#D95319','DisplayName','BD-min');
end 
