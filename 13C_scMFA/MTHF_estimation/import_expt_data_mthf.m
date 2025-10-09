function [MID,SD,time_points,data_metab_list,data_i_list,unlabeled_metabs] = import_expt_data_mthf(input_data_file,min_sd)

[MID,text] = xlsread(input_data_file,'MID','','basic');
[SD,~] = xlsread(input_data_file,'SD','','basic');
[~,unlabeled_metabs] = xlsread(input_data_file,'unlabeled_metabs','','basic');

unlabeled_metabs = string(unlabeled_metabs)';

% save the experimental data in the required format
time_points = text(1,3:end);
data_i_list = MID(2:end,1);
data_metab_list = string(text(2:end,1));
MID = MID(1:end, 2:end);
MID = MID/100;
SD = SD(1:end, 2:end);
SD = SD/100;

% Set a minimum SD value
SD(SD < min_sd) = min_sd;

end