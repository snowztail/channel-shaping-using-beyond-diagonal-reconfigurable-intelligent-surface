clc; clear; close; setup;

[transmit.antenna, reflect.antenna, receive.antenna] = deal(2, 2 .^ (4 : 2 : 6), 2);
[channel.pathloss.direct, channel.pathloss.forward, channel.pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
[number.bond, number.antenna, number.realization, flag.direct] = deal(log2(reflect.antenna) + 1, length(reflect.antenna), 1e2, false);

for r = 1 : number.realization
	% * No RIS
	channel.direct = flag.direct * sqrt(channel.pathloss.direct) * fading_nlos(receive.antenna, transmit.antenna);
	channel.power.direct(r) = norm(channel.direct, 'fro') ^ 2;
	% * Have RIS
	for a = 1 : number.antenna
		channel.forward = sqrt(channel.pathloss.forward) * fading_nlos(reflect.antenna(a), transmit.antenna);
		channel.backward = sqrt(channel.pathloss.backward) * fading_nlos(receive.antenna, reflect.antenna(a));
		for b = 1 : number.bond(a)
			clear scatter_power_max;
			reflect.bond = 2 .^ (b - 1);
			[reflect.beamformer.asymmetric, channel.aggregate.asymmetric] = scatter_power_max(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate.asymmetric(b, a, r) = norm(channel.aggregate.asymmetric, 'fro') ^ 2;

			reflect.beamformer.symmetric.enforced = projection_symmetry(reflect.beamformer.asymmetric);
			channel.aggregate.symmetric.enforced = channel_aggregate(channel.direct, channel.forward, channel.backward, reflect.beamformer.symmetric.enforced);
			channel.power.aggregate.symmetric.enforced(b, a, r) = norm(channel.aggregate.symmetric.enforced, 'fro') ^ 2;

			[reflect.beamformer.symmetric.legacy, channel.aggregate.symmetric.legacy] = scatter_power_max_symmetric_legacy(channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate.symmetric.legacy(b, a, r) = norm(channel.aggregate.symmetric.legacy, 'fro') ^ 2;

			[reflect.beamformer.symmetric.takagi, channel.aggregate.symmetric.takagi] = scatter_power_max_symmetric_takagi(channel.direct, channel.forward, channel.backward, reflect.bond);
			channel.power.aggregate.symmetric.takagi(b, a, r) = norm(channel.aggregate.symmetric.takagi, 'fro') ^ 2;
		end
	end
end
channel.power.aggregate.asymmetric = mean(channel.power.aggregate.asymmetric, ndims(channel.power.aggregate.asymmetric));
channel.power.aggregate.symmetric.enforced = mean(channel.power.aggregate.symmetric.enforced, ndims(channel.power.aggregate.symmetric.enforced));
channel.power.aggregate.symmetric.legacy = mean(channel.power.aggregate.symmetric.legacy, ndims(channel.power.aggregate.symmetric.legacy));
channel.power.aggregate.symmetric.takagi = mean(channel.power.aggregate.symmetric.takagi, ndims(channel.power.aggregate.symmetric.takagi));
save('data/pc_power_symmetry.mat');

figure('Name', 'Channel Power vs RIS Symmetry', 'Position', [0, 0, 500, 400]);
for a = 1 : number.antenna
	handle.power.aggregate(1, a) = plot(1 : number.bond(a), channel.power.aggregate.asymmetric(1 : number.bond(a), a), 'DisplayName', 'Asymmetric: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
	handle.power.aggregate(2, a) = plot(1 : number.bond(a), channel.power.aggregate.symmetric.enforced(1 : number.bond(a), a), 'DisplayName', 'Symmetric-Enforced: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
	hold on;
	handle.power.aggregate(3, a) = plot(1 : number.bond(a), channel.power.aggregate.symmetric.legacy(1 : number.bond(a), a), 'DisplayName', 'Symmetric-Legacy: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
	handle.power.aggregate(4, a) = plot(1 : number.bond(a), channel.power.aggregate.symmetric.takagi(1 : number.bond(a), a), 'DisplayName', 'Symmetric-Takagi: $N_\mathrm{S} = ' + string(reflect.antenna(a)) + '$');
    hold on;
end
hold off;
style_plot(handle.power.aggregate, 4);
legend('Location', 'nw');
hold off; grid on; box on;
set(gca, 'XLim', [1, max(number.bond)], 'XTick', 1 : max(number.bond), 'XTickLabel', '$2^' + string(vec(0 : max(number.bond) - 1)) + '$');
xlabel('RIS Group Size');
ylabel('Channel Power [W]');
savefig('plots/pc_power_symmetry.fig');
matlab2tikz('../assets/simulation/pc_power_symmetry.tex', 'width', '8cm', 'height', '6cm');


function [Theta, H] = scatter_power_max_symmetric_legacy(H_f, H_b, L)
	Theta = eye(size(H_f, 1));
	G = length(Theta) / L;
	for g = 1 : G
		S = (g - 1) * L + 1 : g * L;
        U_f = H_f(S, :) ./ vecnorm(H_f(S, :));
        V_b = H_b(:, S)' ./ vecnorm(H_b(:, S)');
		M = V_b * U_f' + (V_b * U_f').';
		[Q, ~] = takagi_factorization(M);
		Theta(S, S) = Q * Q.';
	end
	H = H_b * Theta * H_f;
end

function [Theta, H] = scatter_power_max_symmetric_takagi(H_d, H_f, H_b, L)
	Theta = eye(size(H_f, 1));

	G = length(Theta) / L;
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	[iter.converge, iter.tolerance, iter.counter, iter.P] = deal(false, 1e-4, 0, norm(H, 'fro') ^ 2);
	while ~iter.converge
		for g = 1 : G
			S = (g - 1) * L + 1 : g * L;
			M = H_b(:, S)' * (H_d + H_b * Theta * H_f) * H_f(S, :)';
			[Q, ~] = takagi_factorization(projection_symmetry(M));
			Theta(S, S) = Q * Q.';
		end
		H = channel_aggregate(H_d, H_f, H_b, Theta);
		P = norm(H, 'fro') ^ 2;
		iter.converge = (abs(P - iter.P) / iter.P <= iter.tolerance);
		iter.P = P;
		iter.counter = iter.counter + 1;
	end
end

function [Q] = projection_symmetry(X)
	Q = (X + X.') / 2;
end
