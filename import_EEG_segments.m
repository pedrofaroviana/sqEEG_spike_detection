%% Importing data segments of equal length from 2-channel EEG
% Requires Fieltrip Toolbox to open edf files
% Will need to run adherence_table.m to extract the edf_files table
% Will return chan X eeg data_segments matrix, and time vector
% Pedro F. Viana, King's College London
% January 2021

% edf_files = adherence_table(pwd,'Europe/London');
% 
% segments_dt = datetime(2019,10,29,10,30,0,'TimeZone','Europe/London'); % example, to replace with annotations vector
% 
% pre_segment_time = 5; % in seconds
% post_segment_time = 5; % in seconds
% 
% % extract srate
% hdr = ft_read_header('/Users/pedrofaroviana/Documents/Neurologia/2019/SUBER/Analysis/dummy_sqEEG.edf');
% srate = hdr.Fs; 

function [data_segments data_segments_time] = import_EEG_segments(path_EEG,edf_files,segments_dt,pre_segment_time,post_segment_time,srate)
ft_warning off
% peri-segment time
peri_segment_time = [pre_segment_time post_segment_time]; % in seconds

for i = 1:size(segments_dt,1)
    disp(['Segment ' num2str(i) '/' num2str(size(segments_dt,1)) ])

A = segments_dt(i) - edf_files.start;
ied_file(i) = find(A == min(A(A>0)));

if i<2 | (i>1 & ied_file(i) ~= ied_file(i-1))
load([path_EEG filesep char(edf_files.filename(ied_file(i)))]);
% ft_read_data(edf_files.filename(ied_file(i)));
time = datenum(edf_files.start(ied_file(i)) + seconds((0:size(data,2)-1)./srate));
else, end

% select only seizure segment
idx = dsearchn(time',datenum(segments_dt(i)));

data_segments(:,:,i) = data(:,idx-round(peri_segment_time(1)*srate):idx+round(peri_segment_time(end)*srate));
data_segments_time(:,:,i) = time(idx-round(peri_segment_time(1)*srate):idx+round(peri_segment_time(end)*srate));

end
ft_warning on
load gong
sound(y,Fs)
end
