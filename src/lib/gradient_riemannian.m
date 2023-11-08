function [G_r] = gradient_riemannian(Theta, G_e)
	G_r = G_e * Theta' - Theta * G_e';
end
