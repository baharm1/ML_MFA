function [MID_metab, MID, SD] = import_experiment_data(file2, min_sd)

%% Read the MID values
[values, MID_metab] = xlsread(file2, 'MID', '', 'basic');
MID = values(:, 2);
SD = values(:, 3);
MID = MID / 100;
SD = SD / 100;
SD(SD < min_sd) = min_sd;




