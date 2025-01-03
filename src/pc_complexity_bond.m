clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(4, [16, 256], 4);
[transmit.power, receive.noise] = deal(db2pow(20), db2pow(-75));
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.antenna, number.realization, flag.direct] = deal(length(reflect.antenna), 1e2, true);

for r = 1 : number.realization
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));

		clear scatter_rate;
		transmit.beamformer.diagonal = precoder_rate(channel.direct, transmit.power, receive.noise);
		[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.diagonal, receive.noise));
		tic;
		while ~iter.converge
			[reflect.beamformer.diagonal, channel.aggregate.diagonal] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer.diagonal, 1, receive.noise);
			transmit.beamformer.diagonal = precoder_rate(channel.aggregate.diagonal, transmit.power, receive.noise);
			receive.rate.diagonal(a, r) = rate_mimo(channel.aggregate.diagonal, transmit.beamformer.diagonal, receive.noise);
			iter.converge = (abs(receive.rate.diagonal(a, r) - iter.rate) / iter.rate <= iter.tolerance);
			iter.rate = receive.rate.diagonal(a, r);
			iter.counter = iter.counter + 1;
		end
		[objective.diagonal(a, r), counter.diagonal(a, r), timer.diagonal(a, r)] = deal(receive.rate.diagonal(a, r), iter.counter, toc);


		clear scatter_rate;
		transmit.beamformer.unitary = precoder_rate(channel.direct, transmit.power, receive.noise);
		[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-4, 0, rate_mimo(channel.direct, transmit.beamformer.unitary, receive.noise));
		tic;
		while ~iter.converge
			[reflect.beamformer.unitary, channel.aggregate.unitary] = scatter_rate(channel.direct, channel.forward, channel.backward, transmit.beamformer.unitary, reflect.antenna(a), receive.noise);
			transmit.beamformer.unitary = precoder_rate(channel.aggregate.unitary, transmit.power, receive.noise);
			receive.rate.unitary(a, r) = rate_mimo(channel.aggregate.unitary, transmit.beamformer.unitary, receive.noise);
			iter.converge = (abs(receive.rate.unitary(a, r) - iter.rate) / iter.rate <= iter.tolerance);
			iter.rate = receive.rate.unitary(a, r);
			iter.counter = iter.counter + 1;
		end
		[objective.unitary(a, r), counter.unitary(a, r), timer.unitary(a, r)] = deal(receive.rate.unitary(a, r), iter.counter, toc);
	end
end

[objective.diagonal, objective.unitary] = deal(mean(objective.diagonal, ndims(objective.diagonal)), mean(objective.unitary, ndims(objective.unitary)));
[counter.diagonal, counter.unitary] = deal(mean(counter.diagonal, ndims(counter.diagonal)), mean(counter.unitary, ndims(counter.unitary)));
[timer.diagonal, timer.unitary] = deal(mean(timer.diagonal, ndims(timer.diagonal)), mean(timer.unitary, ndims(timer.unitary)));

save('data/pc_complexity_bond.mat');
