function [fading] = fading_ricean(component, factor)
	arguments
		component;
		factor = 5;
	end
	fading = sqrt(factor / (1 + factor)) * component.los + sqrt(1 / (1 + factor)) * component.nlos;
end
