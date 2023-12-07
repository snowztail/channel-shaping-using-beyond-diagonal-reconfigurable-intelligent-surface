clc; clear; close; setup;

[transmit.antenna, ris.antenna, receive.antenna] = deal(8, 2 .^ [-inf, 5, 8], 4);
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-75 : -10 : -115));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
transmit.snr = pow2db(transmit.power * channel.pathloss.direct ./ receive.noise);
[number.bond, number.noise, number.antenna, number.realization] = deal(3, length(receive.noise), length(ris.antenna), 1e1);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(ris.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, ris.antenna(a));
		for n = 1 : number.noise
			if ris.antenna(a) == 0
				transmit.beamformer = precoder_rate_pc(channel.direct, transmit.power, receive.noise(n));
				receive.rate(1, n, a, r) = rate_mimo_pc(channel.direct, transmit.beamformer, receive.noise(n));
			else
				ris.bond = [1, 4, ris.antenna(a)];
				for b = 1 : number.bond
					[iter.converge, iter.tolerance, iter.counter] = deal(false, 1e-4, 0);
					ris.scatter = eye(ris.antenna(a));
					channel.aggregate = channel_aggregate(channel.direct, channel.forward, channel.backward, ris.scatter);
					transmit.beamformer = precoder_rate_pc(channel.aggregate, transmit.power, receive.noise(n));
					receive.rate(b, n, a, r) = rate_mimo_pc(channel.aggregate, transmit.beamformer, receive.noise(n));
					while ~iter.converge
						iter.rate = receive.rate(b, n, a, r);
						[ris.scatter, channel.aggregate] = scatter_rate_pc(channel.direct, channel.forward, channel.backward, transmit.beamformer, ris.scatter, ris.bond(b), receive.noise(n));
						transmit.beamformer = precoder_rate_pc(channel.aggregate, transmit.power, receive.noise(n));
						receive.rate(b, n, a, r) = rate_mimo_pc(channel.aggregate, transmit.beamformer, receive.noise(n));
						iter.converge = (abs(receive.rate(b, n, a, r) - iter.rate) / iter.rate <= iter.tolerance);
						iter.counter = iter.counter + 1;
					end
				end
			end
		end
	end
end
receive.rate = mean(receive.rate, 4);
save('data/pc_rate_sx.mat');

figure('Name', 'Achievable Rate vs RIS Configuration', 'Position', [0, 0, 500, 400]);
hold all;
for a = 1 : number.antenna
	if ris.antenna(a) == 0
		handle.rate(1, a) = plot(transmit.snr, receive.rate(1, :, a) / log(2), 'DisplayName', '$N_s = ' + string(ris.antenna(a)) + '$');
	else
		ris.bond = [1, 4, ris.antenna(a)];
		for b = 1 : number.bond
			handle.rate(b, a) = plot(transmit.snr, receive.rate(b, :, a) / log(2), 'DisplayName', '$(N_s, L) = (' + string(ris.antenna(a)) + ', ' + string(ris.bond(b)) + ')$');
		end
	end
end
hold off; grid on;
set(handle.rate(1, find(ris.antenna == 0)), 'Color', 'k', 'Marker', 'none');
style_plot(handle.rate(:, find(ris.antenna ~= 0)), number.bond);
xlabel('Direct SNR [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
legend('Location', 'nw');
savefig('plots/pc_rate_sx.fig');
