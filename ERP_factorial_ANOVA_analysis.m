
%% factorial overall HEP analysis (regardless of the respiratory phase)

% EEGLAB toolbox required
% ERPLAB toolbox required
% MUT toolbox required
% FMUT toolbox required

addpath(genpath('D:\Dati_Andrea\FMUT-master'));
addpath(genpath('D:\Dati_Andrea\massunivariateERPtoolbox'));

% BIN1 = baseline
% BIN2 = respiratory interoception
% BIN3 = cardiac interoception
% BIN4 = respiratory exteroception
% BIN5 = cardiac exteroception

% fronto-central ROI

GND = erplab2GND('gui',...
                 'exclude_chans', {'Pz','P1','P2','P3','P4','CPz','CP1','CP2','CP3','CP4','Cz','C1','C2','C3','C4','F5','F6','FC5','FC6','C5','C6','CP5','CP6','P5','P6','PO3','POz','PO4','O1','O2','Oz','AF7', 'AF8', 'F7', 'F8', 'FT7', 'FT8', 'T7', 'T8', 'TP7', 'TP8', 'P7', 'P8', 'PO7', 'PO8', 'TP9', 'TP10', 'ECG1','ECG2', 'E65', 'E66', 'E67', 'E68', 'E69', 'E70'},...
                 'bsln', [NaN], ...
                 'out_fname', 'no save');

% add bins for following up main effects of condition (interoception vs. exteroception)

GND = bin_mean(GND, 2, 3, 'interoception'); % bin 6
GND = bin_mean(GND, 4, 5, 'exteroception'); % bin 7

% add bins for following up main effects of task (respiratory vs. cardiac)

GND = bin_mean(GND, 2, 4, 'respiratory'); % bin 8
GND = bin_mean(GND, 3, 5, 'cardiac'); % bin 9

% define some variables

n_perm = 10000;
time_wind = [350 597.6];
chan_hood = 0.61;

% full factorial design

GND = FclustGND(GND, ...
                'bins', 2:5, ...  
                'factor_names', {'task', 'condition'}, ...
                'factor_levels', [2, 2], ...
                'time_wind', time_wind, ...
                'mean_wind', 'yes', ...
                'chan_hood', chan_hood, ...
                'n_perm', n_perm, ...
                'save_GND', 'no');

%% factorial HEP analysis (considering respiratory phases)

% BIN1  = baseline inhale
% BIN2  = baseline exhale
% BIN3  = respiratory interoception inhale
% BIN4  = respiratory interoception exhale
% BIN5  = cardiac interoception inhale
% BIN6  = cardiac interoception exhale
% BIN7  = respiratory exteroception inhale
% BIN8  = respiratory exteroception exhale
% BIN9  = cardiac exteroception inhale
% BIN10 = cardiac exteroception exhale

% frontal ROI

GND = erplab2GND('gui',...
                 'include_chans', {'F5','Fp1','Fpz','Fp2','AF3','AFz','AF4','F3','F1','Fz','F2','F4'},...
                 'bsln', [NaN], ...
                 'out_fname', 'no save');

% add bins for following up main effects of condition (interoception vs. exteroception)

GND = bin_mean(GND, 3, 4, 5, 6, 'interoception'); % bin 11
GND = bin_mean(GND, 7, 8, 9, 10, 'exteroception'); % bin 12

% add bins for following up main effects of task (respiratory vs. cardiac)

GND = bin_mean(GND, 3, 4, 7, 8, 'respiratory'); % bin 13
GND = bin_mean(GND, 5, 6, 9, 10, 'cardiac'); % bin 14

% add bins for following up main effects of phase (inhale vs. exhale)

GND = bin_mean(GND, 3, 5, 7, 9, 'inhale'); % bin 15
GND = bin_mean(GND, 4, 6, 8, 10, 'exhale'); % bin 16

% define some variables

n_perm = 10000;
time_wind = [350 597.6];
chan_hood = 0.61;

% full factorial design

GND = FclustGND(GND, ...
                'bins', 3:10, ...  
                'factor_names', {'phase', 'task', 'condition'}, ...  
                'factor_levels', [2, 2, 2], ...
                'time_wind', time_wind, ...
                'mean_wind', 'yes', ...
                'chan_hood', chan_hood, ...
                'n_perm', n_perm, ...
                'save_GND', 'no');

