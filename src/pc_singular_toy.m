clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, reflect.sample] = deal(2, 2, 2, 1e2);
[reflect.azimuth, reflect.elevation] = meshgrid(2 * pi * (1 : reflect.sample) / reflect.sample);

channel.forward = fading_nlos(reflect.antenna, transmit.antenna);
channel.backward = fading_nlos(receive.antenna, reflect.antenna);
for s1 = 1 : reflect.sample
	for s2 = 1 : reflect.sample
		channel.singular.diagonal(s1, s2, :) = svd(channel.backward * scatter_2x2_diagonal(reflect.azimuth(s1, s2), reflect.elevation(s1, s2)) * channel.forward);
		channel.singular.symmtary(s1, s2, :) = svd(channel.backward * scatter_2x2_symmtary(reflect.azimuth(s1, s2), reflect.elevation(s1, s2)) * channel.forward);
	end
end
save('data/pc_singular_toy.mat');

figure('Name', 'Channel Singular Value vs RIS Configuration', 'Position', [0, 0, 500, 400]);
handle.window = tiledlayout(2, 1);
handle.axis(1) = nexttile;
for s = 1 : 2
	handle.mesh(1, s) = mesh(reflect.azimuth, reflect.elevation, channel.singular.diagonal(:, :, s), 'DisplayName', '$\sigma_' + string(s) + '(\mathbf{H})$');
	datatip(handle.mesh(1, s), 'DataIndex', find(vec(channel.singular.diagonal(:, :, s)) == max(vec(channel.singular.diagonal(:, :, s))), 1));
	datatip(handle.mesh(1, s), 'DataIndex', find(vec(channel.singular.diagonal(:, :, s)) == min(vec(channel.singular.diagonal(:, :, s))), 1));
	handle.mesh(1, s).DataTipTemplate.DataTipRows = handle.mesh(1, s).DataTipTemplate.DataTipRows(end);
	handle.mesh(1, s).DataTipTemplate.DataTipRows.Label = '\sigma';
	handle.mesh(1, s).DataTipTemplate.DataTipRows.Format = '%.3f';
	hold on;
end
hold off; legend;
xlabel('$\theta_1$');
ylabel('$\theta_2$');
zlabel({'Singular Value'});
set(gca, 'XTick', 0 : pi / 2 : 2 * pi, 'YTick', 0 : pi / 2 : 2 * pi, 'XTickLabel',{'0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'}, 'YTickLabel',{'0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
title('Diagonal RIS');

handle.axis(2) = nexttile;
for s = 1 : 2
	handle.mesh(2, s) = mesh(reflect.azimuth, reflect.elevation, channel.singular.symmtary(:, :, s), 'DisplayName', '$\sigma_' + string(s) + '(\mathbf{H})$');
	datatip(handle.mesh(2, s), 'DataIndex', find(vec(channel.singular.symmtary(:, :, s)) == max(vec(channel.singular.symmtary(:, :, s))), 1));
	datatip(handle.mesh(2, s), 'DataIndex', find(vec(channel.singular.symmtary(:, :, s)) == min(vec(channel.singular.symmtary(:, :, s))), 1));
	handle.mesh(2, s).DataTipTemplate.DataTipRows = handle.mesh(2, s).DataTipTemplate.DataTipRows(end);
	handle.mesh(2, s).DataTipTemplate.DataTipRows.Label = '\sigma';
	handle.mesh(2, s).DataTipTemplate.DataTipRows.Format = '%.3f';
	hold on;
end
hold off; legend;
xlabel('$\alpha$');
ylabel('$\psi$');
zlabel({'Singular Value'});
set(gca, 'XTick', 0 : pi / 2 : 2 * pi, 'YTick', 0 : pi / 2 : 2 * pi, 'XTickLabel',{'0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'}, 'YTickLabel',{'0', '0.5$\pi$', '$\pi$', '1.5$\pi$', '2$\pi$'});
title('Symmetry Unitary RIS');
linkaxes(handle.axis);
savefig('plots/pc_singular_toy.fig');


function [Theta] = scatter_2x2_diagonal(theta_1, theta_2)
	Theta = diag([exp(1i * theta_1), exp(1i * theta_2)]);
end

function [Theta] = scatter_2x2_symmtary(alpha, psi)
	% * 2x2 unitary symmetry matrix has 3 phase parameters; when direct link is absent, common phase phi has no impact on channel singular value
	phi = 0;
	Theta = exp(1i * phi) * [exp(1i * alpha) * cos(psi), exp(1i * pi / 2) * sin(psi); -exp(-1i * pi / 2) * sin(psi), exp(-1i * alpha) * cos(psi)];
end
