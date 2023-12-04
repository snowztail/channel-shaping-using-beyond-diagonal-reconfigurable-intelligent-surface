clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna] = deal(4, 256, 2);
ris.bond = 2 .^ (0 : 2 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -105));
% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.noise, number.realization] = deal(length(ris.bond), length(receive.noise), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_ricean(receive.antenna, 'ula', transmit.antenna, 'ula');
	% channel.direct = 0 * channel.direct;
	channel.forward = sqrt(channel.pathloss.forward) * fading_ricean(ris.antenna, 'upa', transmit.antenna, 'ula');
	channel.backward = sqrt(channel.pathloss.backward) * fading_ricean(receive.antenna, 'ula', ris.antenna, 'upa');
	channel.sv.direct(:, r) = svd(channel.direct);
	for n = 1 : number.noise
		receive.rate.direct(n, r) = rate_mimo(channel.direct, transmit_max_rate(channel.direct, transmit.power, receive.noise(n)), receive.noise(n));
		for b = 1 : number.bond
			[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-5, 0);
			[ris.scatter, transmit.covariance] = deal(eye(ris.antenna), eye(transmit.antenna) * transmit.power / transmit.antenna);
			channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
			receive.rate.aggregate(b, n, r) = rate_mimo(channel.aggregate, transmit.covariance, receive.noise(n));
			while ~iter.converge
				iter.rate = receive.rate.aggregate(b, n, r);
				[ris.scatter, channel.aggregate] = ris_max_rate(channel.direct, channel.forward, channel.backward, ris.scatter, ris.group(b), transmit.covariance, receive.noise(n));
				transmit.covariance = transmit_max_rate(channel.aggregate, transmit.power, receive.noise(n));
				receive.rate.aggregate(b, n, r) = rate_mimo(channel.aggregate, transmit.covariance, receive.noise(n));
				iter.converge = (abs(receive.rate.aggregate(b, n, r) - iter.rate) <= iter.tolerance);
				iter.counter = iter.counter + 1;
			end
			channel.sv.aggregate(:, b, n, r) = svd(channel.aggregate);
		end
	end
end

channel.sv.direct = mean(channel.sv.direct, 2);
channel.sv.aggregate = mean(channel.sv.aggregate, 4);
receive.rate.direct = mean(receive.rate.direct, 2);
receive.rate.aggregate = mean(receive.rate.aggregate, 3);

figure('Name', 'Achievable Rate vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.rate = gobjects(number.bond + 1, 1);
hold all;
handle.rate(1) = plot(channel.snr, receive.rate.direct / log(2), 'DisplayName', 'No RIS');
for b = 1 : number.bond
	handle.rate(b + 1) = plot(channel.snr, receive.rate.aggregate(b, :) / log(2), 'DisplayName', strcat('$N_g = ', num2str(ris.bond(b)), '$'));
end
hold off; legend('Location', 'nw'); grid on; box on; axis tight;
xlabel('Direct SNR [dB]');
ylabel('Rate [bits/s/Hz]');
style_plot(handle.rate);
savefig('plots/max_rate_bound.fig');


figure('Name', 'Channel Singular Value vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.window = tiledlayout(number.noise, 1);
for n = 1 : number.noise
	handle.axis(n) = nexttile;
	bar([channel.sv.direct, channel.sv.aggregate(:, :, n)]);
	title('Direct SNR $= ' + string(channel.snr(n)) + '$dB');
end
handle.xlabel = xlabel(handle.window, 'Index');
handle.ylabel = ylabel(handle.window, 'Singular Value');
handle.xlabel.Interpreter= 'latex'; handle.ylabel.Interpreter= 'latex';
handle.legend = legend(['No RIS', '$N_g = ' + string(ris.bond) + '$'], 'Orientation', 'horizontal');
handle.legend.ItemTokenSize = [10, 10]; handle.legend.Layout.Tile = 'north';
linkaxes(handle.axis);
savefig('plots/max_rate_bound_sv.fig');
