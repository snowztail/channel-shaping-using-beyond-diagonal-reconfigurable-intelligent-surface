function [W] = transmit_min_interference(H, G)
	W = pagectranspose(receive_min_interference(permute(pagectranspose(H), [1, 2, 4, 3]), pagectranspose(G)));
end
