% set time before and after seizure

function concat_data = concatenate_segments(path_EEG,randomepochs,edf_files,srate)
disp('concatenating data...')
ied_file=zeros(size(randomepochs,1),2);

for i = 1:size(randomepochs,1)
   ied_files(i,1) = min(find(any(isbetween(randomepochs.Var1(i):minutes(1):randomepochs.Var2(i),edf_files.start,edf_files.end),2)));
   ied_files(i,2) = max(find(any(isbetween(randomepochs.Var1(i):minutes(1):randomepochs.Var2(i),edf_files.start,edf_files.end),2)));
end
% 

concat_data = NaN(2,round(3600*srate),size(ied_files,1));

for ri = 1:size(ied_files,1)
    
    time = (edf_files.start(ied_files(ri,1))-hours(1):seconds(1/srate):edf_files.end(ied_files(ri,2))+hours(1))';
    eeg = NaN(2,length(time));

      for fi=1:ied_files(ri,2)-ied_files(ri,1)+1

          idx_datafile = dsearchn(datenum(time),datenum(edf_files.start(ied_files(ri,1)+fi-1)));

    load([path_EEG filesep char(edf_files.filename(ied_files(ri,1)+fi-1))])
    eeg(:,idx_datafile:idx_datafile+size(data,2)-1) = data;
      end

    data2add=eeg(:,isbetween(time,randomepochs.Var1(ri),randomepochs.Var2(ri)));
    concat_data(:,1:length(data2add),ri) = data2add;
end

end
