
%% HEP planned t-test analysis

% overall HEP analysis (regardless of the respiratory phase)

% EEGLAB toolbox required
% ERPLAB toolbox required
% MUT toolbox required

addpath('D:\Dati_Andrea\massunivariateERPtoolbox');

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

GND = bin_dif(GND,3,5,'interoception-exteroception_cardiac');
GND = bin_dif(GND,2,4,'interoception-exteroception_respiratory');
GND = bin_dif(GND,3,2,'interoception_cardiac-interoception_respiratory');
GND = bin_dif(GND,4,3,'cardiac_exteroception-respiratory_exteroception');

GND = clustGND(GND,8,...
               'time_wind',[350 600], ...
               'mean_wind', 'yes',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

%% HEP analysis (considering respiratory phases)

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

GND = bin_dif(GND,6,5,'exhale-inhale_cardiac_interoception');
GND = bin_dif(GND,4,3,'exhale-inhale_respiratory_interoception');
GND = bin_dif(GND,10,9,'exhale-inhale_cardiac_exteroception');
GND = bin_dif(GND,8,7,'exhale-inhale_respiratory_exteroception');

GND = bin_dif(GND,6,10,'exhale_cardiac_interoception-exhale_cardiac_exteroception');
GND = bin_dif(GND,5,9,'inhale_cardiac_interoception-inhale_cardiac_exteroception');
GND = bin_dif(GND,6,4,'exhale_cardiac_interoception-exhale_respiratory_interoception');
GND = bin_dif(GND,5,3,'inhale_cardiac_interoception-inhale_respiratory_interoception');

GND = clustGND(GND,11, ...
               'time_wind',[350 597.6], ...
               'mean_wind', 'yes',...
               'chan_hood',.61, ...
               'thresh_p',.05, ...
               'save_GND','no', ...
               'n_perm', 10000);

