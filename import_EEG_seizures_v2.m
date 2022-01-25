%% Import EEG seizure times
% Read annotations file with seizure times - returns a table (sz_table)
% Pedro F. Viana
% Richardson Lab
% King's College London

% Input should be subject (defined in the folder name), path do txt
% file, start and end dates of study in datetime format
% outputs are diary events for seizures in datetime and datenum formats,
% other events in datetime and datenum formats.
% opts = 'Seizure' for seizures,  'all' for all events



function sz_table = import_EEG_seizures_v2(subject,path_EEG_review,timezone,opts)

% read diary
% cd(path_EEG_review);
file = dir([path_EEG_review filesep subject '_' '*.txt']);
sz_table = readtable([path_EEG_review filesep file.name]);
sz_table = sz_table(:,1:4);

sz_table.Properties.VariableNames = {'type','onset','offset','description'};

is_sz=[];
for i  =  1:size(opts,2)
is_sz(:,i) = startsWith(sz_table.description,opts(i));
end

sz_table = sz_table(logical(sum(is_sz,2)),:);

sz_table.onset = datetime(sz_table.onset,'TimeZone',timezone);
sz_table.offset = datetime(sz_table.offset,'TimeZone',timezone);
sz_table = sortrows(sz_table,'onset');

end


