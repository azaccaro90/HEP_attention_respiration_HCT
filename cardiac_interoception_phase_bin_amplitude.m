
% cardiac interoception phase-bin amplitude

% EEGLAB toolbox required
% CircStat toolbox required

% open all EEG datasets "name_interoception_prunedICA.set"

table_HEP_trials_mean = [];

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

resp = ALLEEG(z).data(66,:);
srate = ALLEEG(z).srate;

% smooth and standardize

resp_elab = filloutliers(resp, 'linear', 'movmedian', srate);
resp_elab_smooth = smoothdata(resp_elab, 'sgolay', srate);
resp_z = zscore(resp_elab_smooth);

% define parameters

% Grund et al., 2022
min_peak_dist = 2; % in sec
frac_iqr = 0.3;   % <--- change this parameter
min_peak_prom = iqr(resp_z) * frac_iqr;

% find peaks and troughts

[pks, pklocs] = findpeaks (resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);
[troughs, trlocs] = findpeaks (-resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);

pktr = [pklocs(:), pks(:); trlocs(:), -troughs(:)];
pktrs = sortrows(pktr);

% always begin with a trought (inhalation)

while pktrs(1,2)>0
      pktrs(1,:) = [];
end

% always end with a peak (exhalation)

while pktrs(end,2)<0
      pktrs(end,:) = [];
end

% always alternate peaks and troughts

for i = 1:4

if pktrs(end,2)>0 && pktrs(end-1,2)>0
         pktrs(end-1,:) = [];
elseif pktrs(end,2)<0 && pktrs(end-1,2)<0
         pktrs(end-1,:) = [];
end

end

i = length(pktrs);

while i > 1
      if pktrs(i,2)>0 && pktrs(i-1,2)>0
         pktrs(i-1,:) = [];
         i = i+1;
      end
      if pktrs(i,2)<0 && pktrs(i-1,2)<0
         pktrs(i-1,:) = [];
         i = i+1;
      end
      i = i-1;
end

% define respiratory cycles

respiration = [];
respiration(:,1) = pktrs([1:2:end],1);                % inhale onset
respiration(:,3) = pktrs([2:2:end],1);                % exhale onset
respiration(:,2) = respiration(:,3)-1;                % inhale offset
respiration([1:end-1],4) = respiration([2:end],1)-1;  % exhale offset
respiration(end,4) = length(resp)-1;
respiration = respiration/srate; % in sec

% delete outlier cycles

respiration(:,5) = respiration(:,4) - respiration(:,1); % cycle duration
list_outliers = isoutlier(respiration);
outliers = find(list_outliers(:,5) == 1);
respiration(outliers,:) = [];
rejected_breaths = length(outliers);
respiration(:,5) = []; % delete cycle duration
inhaleexhale = respiration;

% take data breaks/boundaries into account

event = {ALLEEG(z).event.type}.';
event = str2double(event);
boundaries_list = find(event == -99);
boundaries = ALLEEG(z).event(boundaries_list);
boundaries = [boundaries.latency].';
boundaries = boundaries / srate; % in sec

% delete respiratory cycle if crosses boundaries

i = 1;

while i <= length(respiration)
      for k = 1:length(boundaries)
          if respiration(i,1) <= boundaries(k) && boundaries(k) <= respiration(i,4)
             respiration(i,:) = [];
             i = i+1;
             rejected_breaths = rejected_breaths + 1;
          end
      end
      i = i+1;
end

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

% calculate instantaneous r-peak frequecies

frequency_heartbeats = [];

for i = 1:length(heartbeats) - 1
    
   diff_heartbeats = ((heartbeats(i + 1)- heartbeats(i))); % R-R differencies in sec
   
   frequency_heartbeats = [frequency_heartbeats; (60 / (diff_heartbeats))]; % instantaneous HR
   
end

frequency_heartbeats = [frequency_heartbeats; mean(frequency_heartbeats)];

event_ecg_frequencies = [heartbeats, frequency_heartbeats];

% delete heartbeat if a boundary is between -100 and +600 around a heartbeat

boundary_events = boundaries;

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

sig_chan = [2, 20, 22, 34, 61]; % cardiac interoception exhale-inhale significant channels

HEPs = [];
ECG = [];

for i = 1:length(event_ecg_frequencies)-1
    
    time_window = [(event_ecg_frequencies(i,1) + 0.35) * ALLEEG(z).srate : (event_ecg_frequencies(i,1) + 0.6) * ALLEEG(z).srate]; % in sample points
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
    end

    end

end

% rename r-peaks events based on respiratory phases (inhale vs. exhale)

for i = 1:length(event_ecg_frequencies)

   for k = 1:length(inhaleexhale)
        
   if inhaleexhale(k,1) < event_ecg_frequencies(i, 1) && event_ecg_frequencies(i, 1) < inhaleexhale(k, 2)
      event_ecg_frequencies_tab(i,7) = {['inhale']};
   end
   
   if inhaleexhale(k,3) < event_ecg_frequencies(i, 1) && event_ecg_frequencies(i, 1) < inhaleexhale(k, 4)
      event_ecg_frequencies_tab(i,7) = {['exhale']};
   end
   
   end

end

% delete emply rows

row_has_empty = any(cellfun(@isempty, event_ecg_frequencies_tab(:,[1:7])), 2);
event_ecg_frequencies_tab(row_has_empty,:) = []; % delete
event_ecg_frequencies(row_has_empty,:) = []; % delete

% calculate angular values for each heartbeat

angular_values = [];

for i = 1:length(event_ecg_frequencies)
    
    subsequent = find(inhaleexhale(:,1) > event_ecg_frequencies(i,1));
    previous = find(inhaleexhale(:,1) < event_ecg_frequencies(i,1));
    
    if   length(subsequent) == 0
         inhale_subsequent = inhaleexhale(end);
    else subsequent = subsequent(1);
         previous = previous(end);
         inhale_subsequent = inhaleexhale(subsequent,1);
         inhale_previous = inhaleexhale(previous,1);
    end

    angular_values = [angular_values; ((event_ecg_frequencies(i,1) - inhale_previous) / (inhale_subsequent - inhale_previous)) * 360];

end

event_ecg_frequencies = [event_ecg_frequencies, angular_values];

% calculate mean HEP for each 60 degrees interval

temp_position = find(event_ecg_frequencies(:,5) >= 0 & event_ecg_frequencies(:,5) <= 60);
HEP_trials_60 = event_ecg_frequencies(temp_position,:);
HEP_trials_60_mean(:,3) = mean(HEP_trials_60(:,3));
HEP_trials_60_mean(:,2) = 60;
HEP_trials_60_mean(:,1) = z;

temp_position = find(event_ecg_frequencies(:,5) > 60 & event_ecg_frequencies(:,5) <= 120);
HEP_trials_120 = event_ecg_frequencies(temp_position,:);
HEP_trials_120_mean(:,3) = mean(HEP_trials_120(:,3));
HEP_trials_120_mean(:,2) = 120;
HEP_trials_120_mean(:,1) = z;

temp_position = find(event_ecg_frequencies(:,5) > 120 & event_ecg_frequencies(:,5) <= 180);
HEP_trials_180 = event_ecg_frequencies(temp_position,:);
HEP_trials_180_mean(:,3) = mean(HEP_trials_180(:,3));
HEP_trials_180_mean(:,2) = 180;
HEP_trials_180_mean(:,1) = z;

temp_position = find(event_ecg_frequencies(:,5) > 180 & event_ecg_frequencies(:,5) <= 240);
HEP_trials_240 = event_ecg_frequencies(temp_position,:);
HEP_trials_240_mean(:,3) = mean(HEP_trials_240(:,3));
HEP_trials_240_mean(:,2) = 240;
HEP_trials_240_mean(:,1) = z;

temp_position = find(event_ecg_frequencies(:,5) > 240 & event_ecg_frequencies(:,5) <= 300);
HEP_trials_300 = event_ecg_frequencies(temp_position,:);
HEP_trials_300_mean(:,3) = mean(HEP_trials_300(:,3));
HEP_trials_300_mean(:,2) = 300;
HEP_trials_300_mean(:,1) = z;

temp_position = find(event_ecg_frequencies(:,5) > 300 & event_ecg_frequencies(:,5) <= 360);
HEP_trials_360 = event_ecg_frequencies(temp_position,:);
HEP_trials_360_mean(:,3) = mean(HEP_trials_360(:,3));
HEP_trials_360_mean(:,2) = 360;
HEP_trials_360_mean(:,1) = z;

HEP_trials_mean = [HEP_trials_60_mean; ...
                   HEP_trials_120_mean; ...
                   HEP_trials_180_mean; ...
                   HEP_trials_240_mean; ...
                   HEP_trials_300_mean; ...
                   HEP_trials_360_mean];

% concatenate

table_HEP_trials_mean = [table_HEP_trials_mean; HEP_trials_mean];

end

% export for statistics

xlswrite('table_HEP_trials_mean_60_degrees.xlsx', table_HEP_trials_mean);

