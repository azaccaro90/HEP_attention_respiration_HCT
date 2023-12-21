
% single trial analysis among tasks overall HEP (regardless of the respiratory phase)
% HCT vs. C-TCT

% EEGLAB toolbox required

% open all EEG datasets "name_interoception_prunedICA.set"
% open all EEG datasets "name_exteroception_prunedICA.set"

table_single_trials = [];

for z = 1:length(ALLEEG) % loop among participants
    
% rename events in numeric format

for i = 1:length(ALLEEG(z).event)

ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 16', '16');      % respiratory task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 32', '32');      % respiratory task offset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 48', '48');      % cardiac task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 64', '64');      % cardiac task ofset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'boundary', '-99'); % data breaks
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'ECG', '88');       % r-peaks

end

[~,index] = sortrows([ALLEEG(z).event.latency].'); ALLEEG(z).event = ALLEEG(z).event(index); clear index

% respiratory analysis

srate = ALLEEG(z).srate;

% take data breaks/boundaries into account

event = {ALLEEG(z).event.type}.';
event = str2double(event);
boundaries_list = find(event == -99);
boundaries = ALLEEG(z).event(boundaries_list);
boundaries = [boundaries.latency].';
boundaries = boundaries / srate; % in sec

% define cardiac and respiratory task condition range
                                
list_event_cell = {ALLEEG(z).event.type}.';
list_event = str2double(list_event_cell);

position_16 = find (list_event == 16);
latency_16 = ALLEEG(z).event(position_16);
latency_16 = [latency_16.latency].';
latency_16 = latency_16 / ALLEEG(z).srate;

position_32 = find (list_event == 32);
latency_32 = ALLEEG(z).event(position_32);
latency_32 = [latency_32.latency].';
latency_32 = latency_32 / ALLEEG(z).srate;

position_64 = find (list_event == 64);
latency_64 = ALLEEG(z).event(position_64);
latency_64 = [latency_64.latency].';
latency_64 = latency_64 / ALLEEG(z).srate;

position_48 = find (list_event == 48);
latency_48 = ALLEEG(z).event(position_48);
latency_48 = [latency_48.latency].';
latency_48 = latency_48 / ALLEEG(z).srate;

respiratory_intervals = [latency_16, latency_32];
cardiac_intervals = [latency_64, latency_48];

position_88 = find (list_event == 88);
latency_88 = ALLEEG(z).event(position_88);
latency_88 = [latency_88.latency].';
latency_88 = latency_88 / ALLEEG(z).srate;
heartbeats = latency_88;

% calculate instantaneous R-peak frequecies

frequency_heartbeats = [];

for i = 1:length(heartbeats) - 1
    
   diff_heartbeats = ((heartbeats(i + 1)- heartbeats(i))); % R-R differencies in sec
   
   frequency_heartbeats = [frequency_heartbeats; (60 / (diff_heartbeats))]; % instantaneous HR
   
end

frequency_heartbeats = [frequency_heartbeats; mean(frequency_heartbeats)];

event_ecg_frequencies = [heartbeats, frequency_heartbeats];

boundary_events = boundaries;

% delete heartbeat if a boundary is between -100 and +600 around a heartbeat

i = 1;

while i < length(event_ecg_frequencies)
    
    for k = 1:length(boundary_events)
        
        if event_ecg_frequencies(i,1) - 0.1 <= boundary_events(k) && boundary_events(k) <= event_ecg_frequencies(i,1) + 0.6
           event_ecg_frequencies(i,:) = [];
        end
        
    end
    i = i+1;
end

% calculate mean HEP

sig_chan = [4, 17, 20, 21, 22, 33, 34, 35, 40, 41, 42, 48, 71]; % significant channels HCT vs. C-TCT

HEPs = [];
ECG = [];

for i = 1:length(event_ecg_frequencies)-1
    
    time_window = [(event_ecg_frequencies(i,1) + 0.352) * ALLEEG(z).srate : (event_ecg_frequencies(i,1) + 0.539) * ALLEEG(z).srate]; % in sample points
    time_window = fix(time_window);
    data_window = ALLEEG(z).data(:,time_window);
    HEP_data = data_window([sig_chan] , :); % significant channels
    HEP = mean(mean(HEP_data));
    E65_data = data_window([65], :); % ECG channel
    E65 = mean(E65_data);

    HEPs = [HEPs; HEP];
    ECG = [ECG; E65];

end

% concatenate

event_ecg_frequencies(end,:) = [];
event_ecg_frequencies = [event_ecg_frequencies, HEPs, ECG];

% delete outliers

list_outliers = isoutlier(event_ecg_frequencies);

outliers_second_column = find(list_outliers(:,2) == 1);
outliers_third_column = find(list_outliers(:,3) == 1);
outliers_fourth_column = find(list_outliers(:,4) == 1);

all_outliers = [outliers_second_column; outliers_third_column; outliers_fourth_column];
all_outliers = sort(all_outliers);

event_ecg_frequencies(all_outliers,:) = [];

% create cell matrix

event_ecg_frequencies_tab = num2cell(event_ecg_frequencies);

% rename r-peaks events based on task condition

for j = 1:length(event_ecg_frequencies)
    
    for k = 1:length(cardiac_intervals)
        
    if  cardiac_intervals(k,1) < event_ecg_frequencies(j,1) && event_ecg_frequencies(j,1) < cardiac_intervals(k,2)
        event_ecg_frequencies_tab(j,5) = {['cardiac']};
        event_ecg_frequencies_tab(j,6) = {['interoception']};
%         event_ecg_frequencies_tab(j,6) = {['exteroception']};
    end

    end
    
end

% delete emply rows

row_has_empty = any(cellfun(@isempty, event_ecg_frequencies_tab(:,[1:6])), 2);
event_ecg_frequencies_tab(row_has_empty,:) = []; % delete

% last colomn becomes ID

for i = 1:length(event_ecg_frequencies_tab)
    
    event_ecg_frequencies_tab {i,7} = z;

end

% concatenate

table_single_trials = [table_single_trials; event_ecg_frequencies_tab];

end

table_single_trials(:,1) = table_single_trials(:,7);
table_single_trials(:,7) = [];

% export tables for statistics

% xlswrite('table_single_trials_cardiac_interoception.xlsx', table_single_trials);
% xlswrite('table_single_trials_cardiac_exteroception.xlsx', table_single_trials);
