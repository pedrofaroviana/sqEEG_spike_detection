%% Spike Data Correlation
% Run template on EEG data, sample-by-sample
% inputs: template (chan X time); minimum rho value; minimum and maximum
% amplitude thresholds; table with edf file information; path to EEG data;
% subject ID, sampling rate.

function spikes_table = spike_data_corr(templates,minimum_rho,minimum_p2p,maximum_p2p,edf_files,path_EEG,subject,srate)
tic
ft_warning off
t= toc;
n_templates = size(templates,1);
% loop over edf files
fi=1; n_spikes=0; covered_data=0;
% while (n_spikes < 200 || covered_data < 31) && fi <= size(edf_files,1)
for fi = 1:size(edf_files,1)
    tic
    disp(['File n. ' num2str(fi) '/' num2str(size(edf_files,1)) ', estimated time remaining: ' num2str(mean(t)*(size(edf_files,1)-fi)/60) ' min.'])
    %     disp([num2str(size(spikes_table,1)) ' spikes so far..'])
    % read mat files 
    load([path_EEG filesep char(edf_files.filename(fi))]);

        % initialize temp table, rho value and peak_2_peak (amplitude) matrices
        spikes_table_t_chan = cell(size(templates));

        for ti=1:n_templates
        rho = zeros(length(data)-length(templates{ti,1})+1,2);
        peak_2_peak = zeros(length(data)-length(templates{ti,1})+1,2);
        
        % loop over channels
        for chi = 1:2
            % matrix with sliding window segments
            
            data_corr = buffer(data(chi,:)',length(templates{ti,1}),length(templates{ti,1})-1,'nodelay');
            % correlation coef.
            rho(:,chi) = corr(data_corr,(templates{ti,chi})');
            
            % amplitude
            peak_2_peak(:,chi) = max(data_corr) - min(data_corr);
           
             % indices for segments with > mininum rho and amp. thresholds
            high_corr_idx = find(rho(:,chi)>minimum_rho & peak_2_peak(:,chi) > minimum_p2p(ti,chi) & peak_2_peak(:,chi) < maximum_p2p(ti,chi));
            spikes_table_t_chan{ti,chi} = table(edf_files.start(fi) + seconds(high_corr_idx./srate), rho(high_corr_idx,chi), (zeros(length(high_corr_idx),1) + chi), peak_2_peak(high_corr_idx,chi));

        end

spikes_table_t{ti,fi} = cat(1,spikes_table_t_chan{ti,:});
spikes_table_t{ti,fi}.Properties.VariableNames = {'dt','rho','chan','p2p'};
spikes_table_t{ti,fi} = remove_overlap(spikes_table_t{ti,fi},templates{ti,1},srate);

ft_warning on

t(fi)=toc;
        end
        n_spikes = size(cat(1,spikes_table_t{ti,:}),1);
        disp(['N. spikes: ' num2str(n_spikes)])
        covered_data = caldays(between(edf_files.start(1),edf_files.start(fi),'days'));
        disp(['Days covered: ' num2str(covered_data)])
        fi=fi+1;
end

load handel
sound(y,Fs)
spikes_table{ti} = cat(1,spikes_table_t{ti,:});

toc
end