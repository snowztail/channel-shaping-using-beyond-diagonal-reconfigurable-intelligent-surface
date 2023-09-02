function [steering] = steering_ula(nAnts, angle)
	arguments
		nAnts;
		angle.azimuth = 2 * pi * rand;
	end
	coordinate.x = (0 : nAnts - 1)';
	steering = exp(-1i * pi * coordinate.x * cos(angle.azimuth));
end
