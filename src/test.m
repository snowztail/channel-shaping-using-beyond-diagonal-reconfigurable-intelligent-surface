clear; close; setup;

[base.antenna, ris.antenna, user.antenna] = deal(8, 256, 4);
[power.transmit, power.noise] = deal(db2pow(-20), db2pow(-90));

% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
[pathloss.direct, pathloss.forward, pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
channel.direct = sqrt(pathloss.direct) * fading_ricean(base.antenna, 'ula', user.antenna, 'ula', 10);
channel.forward = sqrt(pathloss.forward) * fading_ricean(base.antenna, 'ula', ris.antenna, 'upa');
channel.backward = sqrt(pathloss.backward) * fading_ricean(ris.antenna, 'upa', user.antenna, 'ula', 10);

ris.connect = num2cell(2 .^ (0 : log2(ris.antenna)));
ris.scatter = repmat({eye(ris.antenna)}, size(ris.connect));

benchmark.rate = rate_mimo(channel.direct, update_bs(channel.direct, power.transmit, power.noise), power.noise);
benchmark.eigenvalue = eig(channel.direct * channel.direct');

nVariables = length(ris.connect);
rate = zeros(1, nVariables);
eigenvalue = zeros(min(base.antenna, user.antenna), nVariables);

for iVariable = 1 : nVariables
	base.covariance{iVariable} = (power.transmit / base.antenna) * eye(base.antenna);
	[iter.converge, iter.tolerance, iter.counter, iter.rate] = deal(false, 1e-5, 0, 0);
	while ~iter.converge
		[ris.scatter{iVariable}, channel.aggregate{iVariable}] = update_ris(channel.direct, channel.forward, channel.backward, ris.scatter{iVariable}, ris.connect{iVariable}, base.covariance{iVariable}, power.noise);
		[base.covariance{iVariable}, base.allocate{iVariable}] = update_bs(channel.aggregate{iVariable}, power.transmit, power.noise);
		rate(iVariable) = rate_mimo(channel.aggregate{iVariable}, base.covariance{iVariable}, power.noise);
		
		iter.converge = (abs(rate(iVariable) - iter.rate) <= iter.tolerance);
		iter.counter = iter.counter + 1;
		iter.rate = rate(iVariable);
	end
	eigenvalue(:, iVariable) = eig(channel.aggregate{iVariable} * channel.aggregate{iVariable}');
end
