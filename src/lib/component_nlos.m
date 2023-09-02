function [component] = component_nlos(nTxs, nRxs)
	component = sqrt(0.5) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
end
