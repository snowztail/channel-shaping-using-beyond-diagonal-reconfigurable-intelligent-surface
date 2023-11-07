function [G_r] = gradient_riemannian(H_d, H_f, H_b, Theta, Q, P_n)
	H = channel_aggregate(H_d, H_f, H_b, Theta);
	
	% * Gradient on the Euclidean space
	G_e = (H_b') / (eye(size(H, 1)) + H * (Q / P_n) * H') * (H * (Q / P_n) * H_f');
	
	% * Gradient on the Riemannian space
	G_r = G_e * Theta' - Theta * G_e';
end
