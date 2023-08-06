
% cardiac interoceptive accuracy

% EEGLAB toolbox required
% Pan-Tompkin algorithm required

% open all EEG datasets "name_interoception_elab.set"

addpath('D:\Dati_Andrea\pan_tompkin');

interoceptive_accuracy_trial = [];

% r-peaks detection

for z = 1:length(ALLEEG)
                   
ecg = ALLEEG(z).data(65,:);
ecg = double(ecg);
ecg = downsample(ecg,2);
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg,256,0);
R_peaks_time = qrs_i_raw / 256;
R_peaks_time = R_peaks_time';

% rename events in numeric format for ERPLAB

for i = 1:length(ALLEEG(z).event)

ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 16', '16');      % respiratory task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 32', '32');      % respiratory task offset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 48', '48');      % cardiac task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 64', '64');      % cardiac task ofset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'boundary', '-99'); % data breaks
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'ECG', '88');       % r-peaks

end

% define cardiac task condition range

list_event_cell = {ALLEEG(z).event.type}.';
list_event = str2double(list_event_cell);

position_64 = find (list_event == 64);
latency_64 = ALLEEG(z).event(position_64);
latency_64 = [latency_64.latency].';
latency_64 = latency_64 / ALLEEG(z).srate;
begin_cardiac = latency_64;

position_48 = find (list_event == 48);
latency_48 = ALLEEG(z).event(position_48);
latency_48 = [latency_48.latency].';
latency_48 = latency_48 / ALLEEG(z).srate;
end_cardiac = latency_48;

cardiac_condition = [begin_cardiac, end_cardiac];

% associate each R-peak with the respective HCT trial

list_R_block1 = [];
list_R_block2 = [];
list_R_block3 = [];
list_R_block4 = [];
list_R_block5 = [];
list_R_block6 = [];
list_R_block7 = [];
list_R_block8 = [];
list_R_block9 = [];
list_R_block10 = [];
list_R_block11 = [];
list_R_block12 = [];

for i = 1:length(R_peaks_time)
    
    if R_peaks_time(i) > cardiac_condition(1,1) && R_peaks_time(i) < cardiac_condition(1,2)
       list_R_block1 = [list_R_block1; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(2,1) && R_peaks_time(i) < cardiac_condition(2,2)
       list_R_block2 = [list_R_block2; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(3,1) && R_peaks_time(i) < cardiac_condition(3,2)
       list_R_block3 = [list_R_block3; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(4,1) && R_peaks_time(i) < cardiac_condition(4,2)
       list_R_block4 = [list_R_block4; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(5,1) && R_peaks_time(i) < cardiac_condition(5,2)
       list_R_block5 = [list_R_block5; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(6,1) && R_peaks_time(i) < cardiac_condition(6,2)
       list_R_block6 = [list_R_block6; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(7,1) && R_peaks_time(i) < cardiac_condition(7,2)
       list_R_block7 = [list_R_block7; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(8,1) && R_peaks_time(i) < cardiac_condition(8,2)
       list_R_block8 = [list_R_block8; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(9,1) && R_peaks_time(i) < cardiac_condition(9,2)
       list_R_block9 = [list_R_block9; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(10,1) && R_peaks_time(i) < cardiac_condition(10,2)
       list_R_block10 = [list_R_block10; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(11,1) && R_peaks_time(i) < cardiac_condition(11,2)
       list_R_block11 = [list_R_block11; R_peaks_time(i)];
    elseif R_peaks_time(i) > cardiac_condition(12,1) && R_peaks_time(i) < cardiac_condition(12,2)
       list_R_block12 = [list_R_block12; R_peaks_time(i)];
    end
    
end

% calculate the real number of heartbeats for each HCT trial

real_heartbeats_block1 = length(list_R_block1);
real_heartbeats_block2 = length(list_R_block2);
real_heartbeats_block3 = length(list_R_block3);
real_heartbeats_block4 = length(list_R_block4);
real_heartbeats_block5 = length(list_R_block5);
real_heartbeats_block6 = length(list_R_block6);
real_heartbeats_block7 = length(list_R_block7);
real_heartbeats_block8 = length(list_R_block8);
real_heartbeats_block9 = length(list_R_block9);
real_heartbeats_block10 = length(list_R_block10);
real_heartbeats_block11 = length(list_R_block11);
real_heartbeats_block12 = length(list_R_block12);

real_heartbeats = [real_heartbeats_block1; real_heartbeats_block2;...
    real_heartbeats_block3; real_heartbeats_block4; real_heartbeats_block5;
    real_heartbeats_block6; real_heartbeats_block7; real_heartbeats_block8;...
    real_heartbeats_block9; real_heartbeats_block10; real_heartbeats_block11;...
    real_heartbeats_block12];

% compare real with counted heartbeats

counted_heartbeats = readtable('counted_heartbeats.xlsx');
counted_heartbeats = table2array(counted_heartbeats);

heartbeats_comparison = [real_heartbeats, counted_heartbeats(:,z)];

% calculate accuracy

interoceptive_accuracy_block1 = 1 - (abs((heartbeats_comparison(1,1) - heartbeats_comparison(1,2))) / heartbeats_comparison(1,1));
interoceptive_accuracy_block2 = 1 - (abs((heartbeats_comparison(2,1) - heartbeats_comparison(2,2))) / heartbeats_comparison(2,1));
interoceptive_accuracy_block3 = 1 - (abs((heartbeats_comparison(3,1) - heartbeats_comparison(3,2))) / heartbeats_comparison(3,1));
interoceptive_accuracy_block4 = 1 - (abs((heartbeats_comparison(4,1) - heartbeats_comparison(4,2))) / heartbeats_comparison(4,1));
interoceptive_accuracy_block5 = 1 - (abs((heartbeats_comparison(5,1) - heartbeats_comparison(5,2))) / heartbeats_comparison(5,1));
interoceptive_accuracy_block6 = 1 - (abs((heartbeats_comparison(6,1) - heartbeats_comparison(6,2))) / heartbeats_comparison(6,1));
interoceptive_accuracy_block7 = 1 - (abs((heartbeats_comparison(7,1) - heartbeats_comparison(7,2))) / heartbeats_comparison(7,1));
interoceptive_accuracy_block8 = 1 - (abs((heartbeats_comparison(8,1) - heartbeats_comparison(8,2))) / heartbeats_comparison(8,1));
interoceptive_accuracy_block9 = 1 - (abs((heartbeats_comparison(9,1) - heartbeats_comparison(9,2))) / heartbeats_comparison(9,1));
interoceptive_accuracy_block10 = 1 - (abs((heartbeats_comparison(10,1) - heartbeats_comparison(10,2))) / heartbeats_comparison(10,1));
interoceptive_accuracy_block11 = 1 - (abs((heartbeats_comparison(11,1) - heartbeats_comparison(11,2))) / heartbeats_comparison(11,1));
interoceptive_accuracy_block12 = 1 - (abs((heartbeats_comparison(12,1) - heartbeats_comparison(12,2))) / heartbeats_comparison(12,1));

interoceptive_accuracy_blocks = [interoceptive_accuracy_block1; interoceptive_accuracy_block2;...
    interoceptive_accuracy_block3; interoceptive_accuracy_block4; interoceptive_accuracy_block5;
    interoceptive_accuracy_block6; interoceptive_accuracy_block7; interoceptive_accuracy_block8;...
    interoceptive_accuracy_block9; interoceptive_accuracy_block10; interoceptive_accuracy_block11;...
    interoceptive_accuracy_block12];

interoceptive_accuracy_trial = [interoceptive_accuracy_trial, interoceptive_accuracy_blocks];

interoceptive_accuracy(z) = mean(interoceptive_accuracy_blocks);

end

