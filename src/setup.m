addpath(genpath(pwd));

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultaxesfontname', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');
set(groot, 'defaultlinelinewidth', 2);
set(groot, 'defaultstemlinewidth', 2);
set(groot, 'defaultscatterlinewidth', 2);

% PBS_ARRAY_INDEX is environment variable (0 -> local, positive -> HPC instance)
disp(getenv('PBS_ARRAY_INDEX')); % rng(str2double(getenv('PBS_ARRAY_INDEX')));
rng shuffle;
