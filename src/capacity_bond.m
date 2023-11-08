clc; clear; close; setup;

[base.antenna, ris.antenna, user.antenna] = deal(4, 8, 2);
ris.bond = 2 .^ (1 : 2 : log2(ris.antenna));
ris.group = ris.antenna ./ ris.bond;
[base.power, user.noise] = deal(db2pow(-20), db2pow(-65 : -20 : -105));
% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.snr = pow2db(base.power * channel.pathloss.direct ./ user.noise);
[number.bond, number.noise, number.realization] = deal(length(ris.bond), length(user.noise), 1e2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_ricean(base.antenna, 'ula', user.antenna, 'ula');
	channel.forward = sqrt(channel.pathloss.forward) * fading_ricean(base.antenna, 'ula', ris.antenna, 'upa');
	channel.backward = sqrt(channel.pathloss.backward) * fading_ricean(ris.antenna, 'upa', user.antenna, 'ula');
	for n = 1 : number.noise
		% benchmark.rate{n, r} = rate_mimo(channel.direct{r}, update_bs(channel.direct{r}, power.transmit, user.noise(n)), user.noise(n));
		for b = 1 : number.bond
			[ris.scatter, base.covariance] = deal(eye(ris.antenna), eye(base.antenna) * base.power / base.antenna);
			[iter.converge, iter.counter, iter.rate, iter.tolerance] = deal(false, 0, 0, 1e-5);
			while ~iter.converge
				[ris.scatter, channel.aggregate] = ris_capacity(channel.direct, channel.forward, channel.backward, ris.scatter, ris.group(b), base.covariance, user.noise(n));
				[ris.scatter, channel.aggregate] = update_ris(channel.direct, channel.forward, channel.backward, ris.scatter, ris.group(b), base.covariance, user.noise(n));
				[base.covariance, base.allocate{b, n, r}] = update_bs(channel.aggregate, power.transmit, user.noise(n));
				user.rate{b, n, r} = rate_mimo(channel.aggregate, base.covariance, user.noise(n));

				iter.converge = (abs(user.rate{b, n, r} - iter.rate) <= iter.tolerance);
				iter.counter = iter.counter + 1;
				iter.rate = user.rate{b, n, r};
			end
			channel.sv{b, n, r} = svd(channel.aggregate);
		end
	end
end

benchmark.rate = mean(cell2mat(benchmark.rate), 2);
benchmark.sv = mean(cell2mat(benchmark.sv), 2);

result.rate = mean(cell2mat(user.rate), 3);
result.sv = squeeze(mean(cell2mat(cellfun(@(x) shiftdim(x, -3), channel.sv, 'UniformOutput', false)), 3));

figure('Name', 'Achievable Rate vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.rate = gobjects(number.bond + 1, 1);
hold all;
handle.rate(1) = plot(pow2db(power.transmit * pathloss.direct ./ cell2mat(user.noise)), benchmark.rate / log(2), 'DisplayName', 'No RIS');
for c = 1 : number.bond
	handle.rate(c + 1) = plot(pow2db(power.transmit * pathloss.direct ./ cell2mat(user.noise)),  result.rate(c, :) / log(2), 'DisplayName', strcat('$N_g = ', num2str(ris.bond(b)), '$'));
end
hold off; legend('Location', 'nw'); grid on; box on; axis tight;
xlabel('Direct SNR [dB]');
ylabel('Rate [bits/s/Hz]');
style_plot(handle.rate);
savefig('plots/rate_connect_rayleigh.fig');


figure('Name', 'Channel Singular Value vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.window = tiledlayout(number.noise, 1);
for n = 1 : number.noise
	handle.axis(n) = nexttile;
	bar([benchmark.sv, shiftdim(result.sv(:, n, :), 2)]);
	title('Direct SNR = ' + string(pow2db(power.transmit * pathloss.direct / user.noise(n))) + 'dB');
end
handle.xlabel = xlabel(handle.window, 'Index');
handle.ylabel = ylabel(handle.window, 'Singular Value');
handle.xlabel.Interpreter= 'latex';
handle.ylabel.Interpreter= 'latex';
handle.legend = legend(['No RIS'; '$N_g = ' + string(ris.bond) + '$'], 'Orientation', 'Horizontal');
handle.legend.ItemTokenSize = [10, 10]; handle.legend.Layout.Tile = 'north';
linkaxes(handle.axis);
savefig('plots/sv_connect_rayleigh.fig');
