function [Theta, H] = scatter_procrustes(H_d, H_f, H_b, opt)
	if opt == "left"
		[U, ~, V] = svd(pinv(H_b) * H_d * H_f');
	elseif opt == "right"
		[U, ~, V] = svd(H_b' * H_d * pinv(H_f));
	end
	Theta = U * V';
	H = channel_aggregate(H_d, H_f, H_b, Theta);
end
