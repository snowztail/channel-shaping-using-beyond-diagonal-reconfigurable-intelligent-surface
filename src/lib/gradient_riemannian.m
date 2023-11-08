function [G_r] = gradient_riemannian(Theta, G_e, S)
	G_r(S, S) = G_e(S, S) * Theta(S, S)' - Theta(S, S) * G_e(S, S)';
end
