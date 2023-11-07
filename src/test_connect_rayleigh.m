clear; close; setup;

[base.antenna, ris.antenna, user.antenna] = deal(4, 64, 4);
[power.transmit, power.noise] = deal(db2pow(-20), num2cell(db2pow(-65 : -20 : -105)));
ris.connect = num2cell(2 .^ (0 : 2 : log2(ris.antenna))');

% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[pathloss.direct, pathloss.forward, pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));

[number.connect, number.noise, number.realization] = deal(length(ris.connect), length(power.noise), 5);
for r = 1 : number.realization
	channel.direct{r} = sqrt(pathloss.direct) * fading_ricean(base.antenna, 'ula', user.antenna, 'ula');
	channel.forward{r} = sqrt(pathloss.forward) * fading_ricean(base.antenna, 'ula', ris.antenna, 'upa');
	channel.backward{r} = sqrt(pathloss.backward) * fading_ricean(ris.antenna, 'upa', user.antenna, 'ula');
	benchmark.sv{r} = svd(channel.direct{r});
	for n = 1 : number.noise
		benchmark.rate{n, r} = rate_mimo(channel.direct{r}, update_bs(channel.direct{r}, power.transmit, power.noise{n}), power.noise{n});
		for c = 1 : number.connect
			[ris.scatter{c, n, r}, base.covariance{c, n, r}] = deal(eye(ris.antenna), (power.transmit / base.antenna) * eye(base.antenna));
			[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-5, 0, 0);
			while ~iter.converge
				[ris.scatter{c, n, r}, channel.aggregate{c, n, r}] = update_ris(channel.direct{r}, channel.forward{r}, channel.backward{r}, ris.scatter{c, n, r}, ris.connect{c}, base.covariance{c, n, r}, power.noise{n});
				[base.covariance{c, n, r}, base.allocate{c, n, r}] = update_bs(channel.aggregate{c, n, r}, power.transmit, power.noise{n});
				user.rate{c, n, r} = rate_mimo(channel.aggregate{c, n, r}, base.covariance{c, n, r}, power.noise{n});
				
				iter.converge = (abs(user.rate{c, n, r} - iter.rate) <= iter.tolerance);
				iter.counter = iter.counter + 1;
				iter.rate = user.rate{c, n, r};
			end
			channel.sv{c, n, r} = svd(channel.aggregate{c, n, r});
		end
	end
end

benchmark.rate = mean(cell2mat(benchmark.rate), 2);
benchmark.sv = mean(cell2mat(benchmark.sv), 2);

result.rate = mean(cell2mat(user.rate), 3);
result.sv = squeeze(mean(cell2mat(cellfun(@(x) shiftdim(x, -3), channel.sv, 'UniformOutput', false)), 3));

figure('Name', 'Achievable Rate vs RIS Group Size', 'Position', [0, 0, 500, 400]);
handle.rate = gobjects(number.connect + 1, 1);
hold all;
handle.rate(1) = plot(pow2db(power.transmit * pathloss.direct ./ cell2mat(power.noise)), benchmark.rate / log(2), 'DisplayName', 'No RIS');
for c = 1 : number.connect
	handle.rate(c + 1) = plot(pow2db(power.transmit * pathloss.direct ./ cell2mat(power.noise)),  result.rate(c, :) / log(2), 'DisplayName', strcat('$N_g = ', num2str(ris.connect{c}), '$'));
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
	title('Direct SNR = ' + string(pow2db(power.transmit * pathloss.direct / power.noise{n})) + 'dB');
end
handle.xlabel = xlabel(handle.window, 'Index');
handle.ylabel = ylabel(handle.window, 'Singular Value');
handle.xlabel.Interpreter= 'latex';
handle.ylabel.Interpreter= 'latex';
handle.legend = legend(['No RIS'; '$N_g = ' + string(ris.connect) + '$'], 'Orientation', 'Horizontal');
handle.legend.ItemTokenSize = [10, 10]; handle.legend.Layout.Tile = 'north';
linkaxes(handle.axis);
savefig('plots/sv_connect_rayleigh.fig');
