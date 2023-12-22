clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna, transmit.stream, network.pair] = deal(4, 2 .^ [4, 6], 4, 3, 2 : 2 : 10);
[transmit.power, receive.noise] = deal(db2pow(-20), db2pow(-100));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.antenna, number.pair, number.realization] = deal(3, length(reflect.antenna), length(network.pair), 20);

for r = 1 : number.realization
	% * No RIS
	for p = 1 : number.pair
		receive.weight = ones(network.pair(p), 1);
		channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair(p), network.pair(p));
		transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
		transmit.beamformer = precoder_rate_bisection_ic(channel.direct, transmit.beamformer, transmit.power, receive.noise, receive.weight);
		receive.wsr.direct(p, r) = receive.weight' * rate_mimo(channel.direct, transmit.beamformer, receive.noise);
	end
	% * Have RIS
	for a = 1 : number.antenna
		reflect.bond = [1, 4, reflect.antenna(a)];
		for p = 1 : number.pair
			receive.weight = ones(network.pair(p), 1);
			channel.direct = sqrt(channel.pathloss.direct) .* fading_nlos(receive.antenna, transmit.antenna, network.pair(p), network.pair(p));
			channel.forward = sqrt(channel.pathloss.forward) .* fading_nlos(reflect.antenna(a), transmit.antenna, 1, network.pair(p));
			channel.backward = sqrt(channel.pathloss.backward) .* fading_nlos(receive.antenna, reflect.antenna(a), network.pair(p), 1);
			transmit.beamformer = precoder_initial_ic(channel.direct, transmit.stream, transmit.power);
			clear scatter_rate_ic;
			for b = 1 : number.bond
				[iter.converge, iter.tolerance, iter.counter, iter.wsr] = deal(false, 1e-3, 0, receive.weight' * rate_mimo(channel.direct, transmit.beamformer, receive.noise));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate] = scatter_rate_ic(channel.direct, channel.forward, channel.backward, transmit.beamformer, reflect.bond(b), receive.noise, receive.weight);
					transmit.beamformer = precoder_rate_bisection_ic(channel.aggregate, transmit.beamformer, transmit.power, receive.noise, receive.weight);
					receive.wsr.aggregate(b, p, a, r) = receive.weight' * rate_mimo(channel.aggregate, transmit.beamformer, receive.noise);
					iter.converge = (abs(receive.wsr.aggregate(b, p, a, r) - iter.wsr) / iter.wsr <= iter.tolerance);
					iter.wsr = receive.wsr.aggregate(b, p, a, r);
					iter.counter = iter.counter + 1;
				end
			end
		end
	end
end
receive.wsr.direct = mean(receive.wsr.direct, ndims(receive.wsr.direct));
receive.wsr.aggregate = mean(receive.wsr.aggregate, ndims(receive.wsr.aggregate));
save('data/ic_rate_user.mat');

figure('Name', 'Weighted Sum-Rate vs Number of Users', 'Position', [0, 0, 500, 400]);
hold all;
handle.wsr.direct = plot(network.pair, receive.wsr.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', '$N^\mathrm{S} = 0$');
for a = 1 : number.antenna
	reflect.bond = [1, 4, reflect.antenna(a)];
	for b = 1 : number.bond
		handle.wsr.aggregate(b, a) = plot(network.pair, receive.wsr.aggregate(b, :, a) / log(2), 'DisplayName', '$(N^\mathrm{S}, L) = (' + string(reflect.antenna(a)) + ', ' + string(reflect.bond(b)) + ')$');
	end
end
style_plot(handle.wsr.aggregate, number.bond);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('User Pairs');
ylabel('Weighted Sum-Rate [bit/s/Hz]');
savefig('plots/ic_rate_user.fig');
matlab2tikz('../assets/simulation/ic_rate_user.tex', 'width', '10cm', 'height', '7.5cm');
