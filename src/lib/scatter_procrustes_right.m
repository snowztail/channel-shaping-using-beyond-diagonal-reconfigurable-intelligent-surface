function [Theta, H] = scatter_procrustes_right(H_d, H_f, H_b)
	[U, ~, V] = svd(H_b' * H_d * pinv(H_f));
	Theta = U * V';
	H = channel_aggregate(H_d, H_f, H_b, Theta);
end
