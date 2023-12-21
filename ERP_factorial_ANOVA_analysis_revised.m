
% HEP planned t-test analysis

% overall HEP analysis regardless of the respiratory phase

addpath('D:\Dati_Andrea\massunivariateERPtoolbox');

% BIN1 = respiratory interoception
% BIN2 = cardiac interoception
% BIN3 = respiratory exteroception
% BIN4 = cardiac exteroception

% exclude external channels

GND = erplab2GND('gui',...
                 'exclude_chans', {'O1','O2','Oz','F7', 'F8', 'FT7', 'FT8', 'T7', 'T8', 'TP7', 'TP8', 'P7', 'P8', 'PO7', 'PO8', 'TP9', 'TP10', 'ECG1','ECG2', 'E65', 'E66', 'E67', 'E68', 'E69', 'E70'},...
                 'bsln', [NaN], ...
                 'out_fname', 'no save');

% downsample from 512hz to 128hz using boxcar filter
GND = decimateGND(GND, 4, 'boxcar', [NaN], 'no', 0);

GND = bin_dif(GND,2,4,'interoception-exteroception_cardiac'); % bin 5
GND = bin_dif(GND,1,3,'interoception-exteroception_respiratory'); % bin 6

% GND = bin_dif(GND,2,1,'interoception_cardiac-interoception_respiratory');
% GND = bin_dif(GND,4,3,'cardiac_exteroception-respiratory_exteroception');

GND = bin_dif(GND,5,6,'task*condition_interaction'); % creates bin7 for interaction

GND = clustGND(GND,7,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

% main effects of condition (interoception vs. exteroception)
             
GND = bin_mean(GND, 1, 2, 'interoception'); % bin 8
GND = bin_mean(GND, 3, 4, 'exteroception'); % bin 9

% main effects of task (respiratory vs. cardiac)

GND = bin_mean(GND, 1, 3, 'respiratory'); % bin 10
GND = bin_mean(GND, 2, 4, 'cardiac'); % bin 11

% create bin for main effects

GND = bin_dif(GND,8,9,'interoception-exteroception'); % bin 12
GND = bin_dif(GND,11,10,'cardiac-respiratory'); % bin 13

% test main effects

GND = clustGND(GND,13,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

%% HEP analysis taking respiratory phases into account

% BIN1  = respiratory interoception inhale
% BIN2  = respiratory interoception exhale
% BIN3  = cardiac interoception inhale
% BIN4  = cardiac interoception exhale
% BIN5  = respiratory exteroception inhale
% BIN6  = respiratory exteroception exhale
% BIN7  = cardiac exteroception inhale
% BIN8 = cardiac exteroception exhale

% fronto-central ROI
           
GND = erplab2GND('gui',...
                 'include_chans', {'AF7','AF8','Fp1','Fp2','Fpz','AF3','F1','FC3','FC1','AFz','Fz','FCz','AF4','F5','F6','F3','F4','F2','FC4','FC2'},...
                 'bsln', [NaN], ...
                 'out_fname', 'no save');
             
% downsample from 512hz to 128hz using boxcar filter
GND = decimateGND(GND, 4, 'boxcar', [NaN], 'no', 0);

GND = bin_dif(GND,4,3,'deltaHEP_HCT'); % bin 9
GND = bin_dif(GND,2,1,'deltaHEP_BCT'); % bin 10
GND = bin_dif(GND,8,7,'deltaHEP_C-TCT'); % bin 11
GND = bin_dif(GND,6,5,'deltaHEP_B-TCT'); % bin 12

% GND = bin_dif(GND,4,8,'exhale_HCT-C-TCT');
% GND = bin_dif(GND,4,2,'exhale_HCT-BCT');
% GND = bin_dif(GND,3,7,'inhale_HCT-C-TCT');
% GND = bin_dif(GND,3,1,'inhale_HCT-BCT');

GND = bin_dif(GND,9,11,'cardiac_interaction'); % bin 13
GND = bin_dif(GND,10,12,'respiratory_interaction'); % bin 14

GND = bin_dif(GND,13,14,'2x2x2_interaction'); % bin 15 for triple interaction

GND = clustGND(GND,15,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

% main effects of condition (interoception vs. exteroception)
             
GND = bin_mean(GND, 1, 2, 3, 4, 'interoception'); % bin 16
GND = bin_mean(GND, 5, 6, 7, 8, 'exteroception'); % bin 17

% main effects of task (respiratory vs. cardiac)

GND = bin_mean(GND, 1, 2, 5, 6, 'respiratory'); % bin 18
GND = bin_mean(GND, 3, 4, 7, 8, 'cardiac'); % bin 19

% main effects of phase (inhale vs. exhale)

GND = bin_mean(GND, 1, 3, 5, 7, 'inhale'); % bin 20
GND = bin_mean(GND, 2, 4, 6, 8, 'exhale'); % bin 21

% create bin for main effects

GND = bin_dif(GND,16,17,'interoception-exteroception'); % bin 22
GND = bin_dif(GND,19,18,'cardiac-respiratory'); % bin 23
GND = bin_dif(GND,21,20,'exhale-inhale'); % bin 24

GND = clustGND(GND,22,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

% 2x2 interactions

GND = bin_dif(GND,23,22,'task*condition_interaction'); % bin 25 task x condition interaction
GND = bin_dif(GND,23,24,'task*phase_interaction'); % bin 26 task x phase interaction
GND = bin_dif(GND,22,24,'condition*phase_interaction'); % bin 27 condition x phase interaction

GND = clustGND(GND,25,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

           
%% Baseline with respiratory phases

% BIN1  = baseline inhale
% BIN2  = baseline exhale

% exclude external channels

GND = erplab2GND('gui',...
                 'exclude_chans', {'O1','O2','Oz','F7', 'F8', 'FT7', 'FT8', 'T7', 'T8', 'TP7', 'TP8', 'P7', 'P8', 'PO7', 'PO8', 'TP9', 'TP10', 'ECG1','ECG2', 'E65', 'E66', 'E67', 'E68', 'E69', 'E70'},...
                 'bsln', [NaN], ...
                 'out_fname', 'no save');

% downsample from 512hz to 128hz using boxcar filter
GND = decimateGND(GND, 4, 'boxcar', [NaN], 'no', 0);

GND = bin_dif(GND,2,1,'deltaHEP_baseline'); % bin 3

GND = clustGND(GND,3,...
               'time_wind',[350 600], ...
               'mean_wind', 'no',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

