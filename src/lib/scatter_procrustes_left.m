function [Theta, H] = scatter_procrustes_left(H_d, H_f, H_b)
	[U, ~, V] = svd(pinv(H_b) * H_d * H_f');
	Theta = U * V';
	H = channel_aggregate(H_d, H_f, H_b, Theta);
end
