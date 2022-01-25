%% Read data headers and extract time start and end
% Fieltrip Toolbox is required to open the edf files
% Pedro F. Viana, King's College London
% January 2021

function edf_files = adherence_table(path_EEG,time_zone)

if strcmp(path_EEG(end-2:end),'mat')

filelist = struct2table(dir([path_EEG filesep '*.mat']));
filelist(startsWith(filelist.name,'._'),:) = [];

edf_files = table(filelist.name,NaT(size(filelist.name),'TimeZone',time_zone),NaT(size(filelist.name),'TimeZone',time_zone),'VariableNames',{'filename','start','end'});

for fi = 1:size(edf_files,1)
    % Read header with Fieldtrip
    load([path_EEG filesep char(filelist.name(fi))],'-mat','hdr')
    edf_files.start(fi) = datetime(hdr.orig.T0,'TimeZone',time_zone);
    edf_files.end(fi) = edf_files.start(fi) + seconds(hdr.nSamples / hdr.Fs);
end

else
ft_warning off
filelist = struct2table(dir([path_EEG filesep '*.edf']));
filelist(startsWith(filelist.name,'._'),:) = [];
edf_files = table(filelist.name,NaT(size(filelist.name),'TimeZone',time_zone),NaT(size(filelist.name),'TimeZone',time_zone),'VariableNames',{'filename','start','end'});

for fi = 1:size(edf_files,1)
    % Read header with Fieldtrip
    hdr = ft_read_header([path_EEG filesep char(edf_files.filename(fi))]);

    % Extract T0 and duration of each file
    edf_files.start(fi) = datetime(hdr.orig.T0,'TimeZone',time_zone);
    edf_files.end(fi) = edf_files.start(fi) + seconds(hdr.nSamples / hdr.Fs);
end
end
% sort by time
edf_files = sortrows(edf_files,2);
ft_warning on

end