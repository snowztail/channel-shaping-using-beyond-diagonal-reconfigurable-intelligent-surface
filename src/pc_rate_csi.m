clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, 128, 4);
[transmit.power, receive.noise] = deal(db2pow(-20 : 5 : 20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[reflect.bond, channel.uncertainty] = deal([1, reflect.antenna], [1e-2, 1e-1, 5e-1]);
[number.bond, number.power, number.uncertainty, number.realization] = deal(length(reflect.bond), length(transmit.power), length(channel.uncertainty), 2);

for r = 1 : number.realization
	channel.direct = sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	% * No RIS
	for p = 1 : number.power
		transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
		receive.rate.direct(p, r) = rate_mimo(channel.direct, transmit.beamformer, receive.noise);
	end
	% * Have RIS
	channel.forward.actual = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna, transmit.antenna);
	channel.backward.actual = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna);
	for u = 1 : number.uncertainty
		channel.forward.estimate = estimation_effect(channel.forward.actual, channel.uncertainty(u));
		channel.backward.estimate = estimation_effect(channel.backward.actual, channel.uncertainty(u));
		clear scatter_rate;
		for b = 1 : number.bond
			for p = 1 : number.power
				transmit.beamformer = precoder_rate(channel.direct, transmit.power(p), receive.noise);
				[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer, receive.noise));
				while ~iter.converge
					[reflect.beamformer, channel.aggregate.estimate] = scatter_rate(channel.direct, channel.forward.estimate, channel.backward.estimate, transmit.beamformer, reflect.bond(b), receive.noise);
					transmit.beamformer = precoder_rate(channel.aggregate.estimate, transmit.power(p), receive.noise);
					receive.rate.estimate(p, b, u, r) = rate_mimo(channel.aggregate.estimate, transmit.beamformer, receive.noise);
					iter.converge = (abs(receive.rate.estimate(p, b, u, r) - iter.rate) / iter.rate <= iter.tolerance);
					iter.rate = receive.rate.estimate(p, b, u, r);
					iter.counter = iter.counter + 1;
				end
				channel.aggregate.actual = channel_aggregate(channel.direct, channel.forward.actual, channel.backward.actual, reflect.beamformer);
				receive.rate.aggregate(p, b, u, r) = rate_mimo(channel.aggregate.actual, transmit.beamformer, receive.noise);
			end
		end
	end
end
receive.rate.direct = mean(receive.rate.direct, ndims(receive.rate.direct));
receive.rate.aggregate = mean(receive.rate.aggregate, ndims(receive.rate.aggregate));
save('data/pc_rate_csi.mat');

figure('Name', 'Achievable Rate vs Channel Estimation Error', 'Position', [0, 0, 500, 400]);
hold all;
handle.rate.direct = plot(pow2db(transmit.power), receive.rate.direct / log(2), 'Color', 'k', 'Marker', 'none', 'DisplayName', 'No RIS');
for u = 1 : number.uncertainty
	handle.rate.aggregate(1, u) = plot(pow2db(transmit.power), receive.rate.aggregate(:, 1, u) / log(2), 'DisplayName', 'D: $\epsilon = ' + string(channel.uncertainty(u)) + '$');
end
for u = 1 : number.uncertainty
	handle.rate.aggregate(2, u) = plot(pow2db(transmit.power), receive.rate.aggregate(:, 2, u) / log(2), 'DisplayName', 'BD: $\epsilon = ' + string(channel.uncertainty(u)) + '$');
end
style_plot(handle.rate.aggregate', number.uncertainty);
hold off; grid on; ylim tight; box on; legend('Location', 'nw');
xlabel('Transmit Power [dB]');
ylabel('Achievable Rate [bit/s/Hz]');
savefig('plots/pc_rate_csi.fig');
matlab2tikz('../assets/simulation/pc_rate_csi.tex', 'width', '10cm', 'height', '7.5cm');



function [H] = estimation_effect(H, sigma)
	H = H .* (1 + sigma * sqrt(0.5) * (randn(size(H)) + 1i * randn(size(H))));
end
