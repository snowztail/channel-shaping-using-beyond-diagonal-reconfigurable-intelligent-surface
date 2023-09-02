function [channel] = create_channel_indoor(nTxs, nRxs, nSxs)
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
