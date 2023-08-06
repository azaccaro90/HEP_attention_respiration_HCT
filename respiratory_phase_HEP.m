%% Respiratory analysis

% EEGLAB toolbox required
% ERPLAB toolbox required

for z = 1:length(ALLEEG) % loop among participants

% import data

resp = ALLEEG(z).data(66,:);
srate = ALLEEG(z).srate;
time = (0:length(resp)-1)/srate;

% smooth and standardize

resp_elab = filloutliers(resp, 'linear', 'movmedian', srate);
resp_elab_smooth = smoothdata(resp_elab, 'sgolay', srate);
resp_z = zscore(resp_elab_smooth);

% check smoothing

% figure(z+1)
% plot(time,resp)
% hold on
% plot(time, resp_elab_smooth)

% define parameters

% Grund et al., 2022
min_peak_dist = 2; % in sec
frac_iqr = 0.3;   % <--- change this parameter
min_peak_prom = iqr(resp_z) * frac_iqr;

% Power et al., 2020
% min_peak_dist = 2; % in sec
% min_peak_prom = 0.5;

% find peaks and troughts

[pks, pklocs] = findpeaks (resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);
[troughs, trlocs] = findpeaks (-resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);

% check peaks

% figure(3)
% findpeaks (resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);
% 
% figure(4)
% findpeaks (-resp_z, 'minpeakdistance', min_peak_dist*srate, 'minpeakprominence', min_peak_prom);

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

% check

for i = 2:length(pktrs)
    if pktrs(i,2)>0 && pktrs(i-1,2)>0
       disp (ALLEEG(z).setname), disp 'has 2 peaks!'
    elseif pktrs(i,2)<0 && pktrs(i-1,2)<0
       disp (ALLEEG(z).setname), disp 'has 2 troughts!'
    end
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

% rename events in numeric format for ERPLAB

for i = 1:length(ALLEEG(z).event)

ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 16', '16');      % respiratory task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 32', '32');      % respiratory task offset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 48', '48');      % cardiac task onset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'S 64', '64');      % cardiac task ofset
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'boundary', '-99'); % data breaks
ALLEEG(z).event(i).type = strrep(ALLEEG(z).event(i).type, 'ECG', '88');       % r-peaks

end

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

% % final check
% 
% respiration_sp = respiration * srate; % in sample points
% inhale_onset  = respiration_sp(:,1);
% inhale_offset = respiration_sp(:,2);
% exhale_onset  = respiration_sp(:,3);
% exhale_offset = respiration_sp(:,4);
% 
% % check inhalation
% 
% figure(5)
% plot(time,resp_z,'-o', 'LineWidth',1,'MarkerIndices', respiration_sp(:,1), 'MarkerSize', 10,  'MarkerEdgeColor','r')
% hold on
% plot(time,resp_z,'-o', 'LineWidth',1,'MarkerIndices', respiration_sp(:,2), 'MarkerSize', 10,  'MarkerEdgeColor','k')
% 
% title('Respiration - Check inhalation')
% xlabel('Time(sec)')
% ylabel('Amplitude(z-score)')
% legend('inhale_onset','inhale_offset')
% 
% % check exhalation
% 
% figure(6)
% plot(time,resp_z,'-o', 'LineWidth',1,'MarkerIndices', respiration_sp(:,3), 'MarkerSize', 10,  'MarkerEdgeColor','r')
% hold on
% plot(time,resp_z,'-o', 'LineWidth',1,'MarkerIndices', respiration_sp(:,4), 'MarkerSize', 10,  'MarkerEdgeColor','k')
% 
% title('Respiration - Check exhalation')
% xlabel('Time(sec)')
% ylabel('Amplitude(z-score)')
% legend('exhale_onset','exhale_offset')

% calculate some respiratory parameters

average_inhale_duration(z) = mean(respiration(:, 2)) - mean(respiration(:, 1));
average_exhale_duration(z) = mean(respiration(:, 4)) - mean(respiration(:, 3));
average_breath_duration(z) = average_inhale_duration(z) + average_exhale_duration(z);
breath_frequency(z) = 60 / average_breath_duration(z);
IE_ratio(z) = average_inhale_duration(z) / average_exhale_duration(z);

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
heartbeats = [position_88, latency_88];

% rename r-peaks events based on task condition

heartbeats_cardiac_condition = []; 
heartbeats_respiratory_condition = [];

for j = 1:length(heartbeats)
    
    for k = 1:length(cardiac_intervals)
        
    if  cardiac_intervals(k,1) < heartbeats(j, 2) && heartbeats(j, 2) < cardiac_intervals(k,2)
        heartbeats_cardiac_condition = [heartbeats_cardiac_condition; heartbeats(j,:)];
    end
    
    end
    
    for k = 1:length(respiratory_intervals)
        
    if  respiratory_intervals(k,1) < heartbeats(j, 2) && heartbeats(j, 2) < respiratory_intervals(k,2)
        heartbeats_respiratory_condition = [heartbeats_respiratory_condition; heartbeats(j,:)];
    end
    
    end
    
end

% rename events

[ALLEEG(z).event(heartbeats_respiratory_condition(:,1)).type] = deal(['88111']); % r-peaks respiratory condition

[ALLEEG(z).event(heartbeats_cardiac_condition(:,1)).type] = deal(['88222']); % r-peaks cardiac condition

% rename r-peaks events based on respiratory phases (inhale vs. exhale)

list_event_cell = {ALLEEG(z).event.type}.';
list_event = str2double(list_event_cell);

heartbeats_respiratory_condition = find (list_event == 88111);
latency_heartbeats_respiratory_condition = ALLEEG(z).event(heartbeats_respiratory_condition);
latency_heartbeats_respiratory_condition = [latency_heartbeats_respiratory_condition.latency].';
latency_heartbeats_respiratory_condition = latency_heartbeats_respiratory_condition / ALLEEG(z).srate;
heartbeats_respiratory_condition = [heartbeats_respiratory_condition, latency_heartbeats_respiratory_condition];

heartbeats_cardiac_condition = find (list_event == 88222);
latency_heartbeats_cardiac_condition = ALLEEG(z).event(heartbeats_cardiac_condition);
latency_heartbeats_cardiac_condition = [latency_heartbeats_cardiac_condition.latency].';
latency_heartbeats_cardiac_condition = latency_heartbeats_cardiac_condition / ALLEEG(z).srate;
heartbeats_cardiac_condition = [heartbeats_cardiac_condition, latency_heartbeats_cardiac_condition];

inhaleexhale = respiration;

% cardiac condition

heartbeats_cardiac_condition_inhale = [];
heartbeats_cardiac_condition_exhale = [];

for i = 1:length(heartbeats_cardiac_condition)

   for k = 1:length(inhaleexhale)
        
   if inhaleexhale(k,1) < heartbeats_cardiac_condition(i, 2) && heartbeats_cardiac_condition(i, 2) < inhaleexhale(k, 2)
      heartbeats_cardiac_condition_inhale = [heartbeats_cardiac_condition_inhale; heartbeats_cardiac_condition(i,:)];
   end
   
   if inhaleexhale(k,3) < heartbeats_cardiac_condition(i, 2) && heartbeats_cardiac_condition(i, 2) < inhaleexhale(k, 4)
      heartbeats_cardiac_condition_exhale = [heartbeats_cardiac_condition_exhale; heartbeats_cardiac_condition(i,:)];
   end
   
   end

end

% respiratory condition

heartbeats_respiratory_condition_inhale = [];
heartbeats_respiratory_condition_exhale = [];

for i = 1:length(heartbeats_respiratory_condition)

   for k = 1:length(inhaleexhale)
        
   if inhaleexhale(k,1) < heartbeats_respiratory_condition(i, 2) && heartbeats_respiratory_condition(i, 2) < inhaleexhale(k, 2)
      heartbeats_respiratory_condition_inhale = [heartbeats_respiratory_condition_inhale; heartbeats_respiratory_condition(i,:)];
   end
   
   if inhaleexhale(k,3) < heartbeats_respiratory_condition(i, 2) && heartbeats_respiratory_condition(i, 2) < inhaleexhale(k, 4)
      heartbeats_respiratory_condition_exhale = [heartbeats_respiratory_condition_exhale; heartbeats_respiratory_condition(i,:)];
   end
   
   end

end

% rename events

[ALLEEG(z).event(heartbeats_respiratory_condition_inhale(:,1)).type] = deal(['881']); % r-peaks respiratory condition inhale
[ALLEEG(z).event(heartbeats_respiratory_condition_exhale(:,1)).type] = deal(['882']); % r-peaks respiratory condition exhale

[ALLEEG(z).event(heartbeats_cardiac_condition_inhale(:,1)).type] = deal(['883']); % r-peaks cardiac condition inhale
[ALLEEG(z).event(heartbeats_cardiac_condition_exhale(:,1)).type] = deal(['884']); % r-peaks cardiac condition exhale

[~,index] = sortrows([ALLEEG(z).event.latency].'); ALLEEG(z).event = ALLEEG(z).event(index); clear index

 end

%% ERPLAB-HEP calculation

% create binlist

for z = 1:length(ALLEEG)

ALLEEG(z) = pop_editeventlist(ALLEEG(z), ...
           'AlphanumericCleaning', 'on', ...
           'BoundaryNumeric', { -99}, ...
           'BoundaryString', { 'boundary' }, ...
           'ExportEL', 'elist_interoception_phases.txt', ...
           'List', 'D:\Dati_Andrea\data\Breath-Heart Interoception_NEW\equation_list_interoception_phases.txt', ...
           'SendEL2', 'All', ...
           'UpdateEEG', 'codelabel', ...
           'Warning', 'off' );

% extract bin-based epochs

ALLEEG(z) = pop_epochbin(ALLEEG(z) , [-100.0  600.0], 'none');

% delete epochs with more than 1 event (double events)

num_event = [];

for i = 1:length(ALLEEG(z).epoch)

    num_event = [num_event; length(ALLEEG(z).epoch(i).eventbini)];

end

double_event = find(num_event > 1);

if ~isempty(double_event)

   ALLEEG(z) = pop_select(ALLEEG(z), 'notrial', double_event); % delete epochs with more than 1 event

end

ALLEEG(z).deleted_epochs = length(double_event);

% correct artifacts

ALLEEG(z) = pop_artmwppth(ALLEEG(z), ...
            'Channel', [ 1:28 31:64 71], ...
            'Flag',  1, ...
            'Threshold',  100, ...
            'Twindow', [-100 600], ...
            'Windowsize', 200, ...
            'Windowstep', 100 );
        
% sync EEGLAB-ERPLAB artifacts

ALLEEG(z) = pop_syncroartifacts(ALLEEG(z), 'Direction', 'eeglab2erplab');

% average epochs and save each .erp files

ALLERP(z) = pop_averager(ALLEEG(z), ...
            'Criterion', 'good', ...
            'DQ_flag', 1, ...
            'ExcludeBoundary', 'on', ...
            'SEM', 'on', ...
            'Saveas', 'on', ...
            'Warning', 'off' );

end

