[nTxs, nSxs, nRxs] = deal(4, 16, 2);
[power.transmit, power.noise] = deal(db2pow(-20), db2pow(-100));
[channel] = create_channel(nTxs, nRxs, nSxs);
[ris] = create_ris(nSxs);


function [channel] = create_channel(nTxs, nRxs, nSxs)
	[pathloss.direct, pathloss.forward, pathloss.backward] = deal(db2pow(-65), db2pow(-54), db2pow(-46));
	% [distance.direct, distance.forward, distance.backward] = deal(-14.7, -10, -6.3);
	% [exponent.direct, exponent.forward, exponent.backward] = deal(-3, -2.4, -2);
	
	[steering.transmit, steering.scatter, steering.receive] = deal(steering_ula(nTxs), steering_upa(nSxs), steering_ula(nRxs));
	
	[fading.direct.los, fading.direct.nlos] = deal(component_los(steering.transmit, steering.receive), component_nlos(nTxs, nRxs));
	[fading.forward.los, fading.forward.nlos] = deal(component_los(steering.transmit, steering.scatter), component_nlos(nTxs, nSxs));
	[fading.backward.los, fading.backward.nlos] = deal(component_los(steering.scatter, steering.receive), component_nlos(nSxs, nRxs));
	
	channel.direct = sqrt(pathloss.direct) * fading_ricean(fading.direct);
	channel.forward = sqrt(pathloss.forward) * fading_ricean(fading.forward);
	channel.backward = sqrt(pathloss.backward) * fading_ricean(fading.backward);
end

function [ris] = create_ris(nSxs)
	% x = randn(nSxs) + 1i * randn(nSxs);
	% ris.scatter = x * (x' * x) ^ -0.5;
	ris.scatter = eye(nSxs);
end

function [steering] = steering_ula(nAnts, angle)
	arguments
		nAnts;
		angle.azimuth = 2 * pi * rand;
	end
	coordinate.x = (0 : nAnts - 1)';
	steering = exp(-1i * pi * coordinate.x * cos(angle.azimuth));
end

function [steering] = steering_upa(nAnts, angle)
	arguments
		nAnts;
		angle.azimuth = 2 * pi * rand;
		angle.elevation = 2 * pi * rand;
	end
	[coordinate.x, coordinate.y] = meshgrid(0 : sqrt(nAnts) - 1);
	coordinate.r = coordinate.x(:) * cos(angle.azimuth) + coordinate.y(:) * sin(angle.azimuth);
	steering = exp(-1i * pi * coordinate.r * cos(angle.elevation));
end

function [component] = component_los(tsv, rsv)
	component = rsv * tsv';
end

function [component] = component_nlos(nTxs, nRxs)
	component = sqrt(0.5) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
end

function [fading] = fading_ricean(component, factor)
	arguments
		component;
		factor = 5;
	end
	fading = sqrt(factor / (1 + factor)) * component.los + sqrt(1 / (1 + factor)) * component.nlos;
end
