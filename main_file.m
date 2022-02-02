% Template Matching Spike Detection with ULT sqEEG
% Pedro Faro Viana, 20/Jan/2022
% IoPPN, King's College London

% subjects to analyse (subject IDs and their timezones)

subjects_IED = {'E02','E04','E06','E08','E09','SUB001','SUB002'; ...
                'Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/London','Europe/London'}

% looping analysis over subjects
for subji = 1:length(subjects_IED)
    close all
subject = subjects_IED{1,subji}

% timezone
time_zone = subjects_IED{2,subji}

% path to EEG annotations file - folder containing .txt files with seizure annotations, extracted from Episight reader
path_EEG_review = '';

% path to EEG data - containing subfolders with each patient data (folder names matching the subject IDs)
path_EEG = '';

% path to ied detection analysis folder
path_IED = '';

% channel labels
channel_label = {'D-C','P-C'};

%% Import adherence and seizure data
% 
% % write a copy of EEG data .edf files to .mat files (in a subfolder) to speed up analysis
% cd([path_EEG filesep subject])
% edf_2_mat

%%

% extract headers, onset and offset times of each EEG file into a table - adherence_table function
edf_files = adherence_table([path_EEG filesep subject filesep 'mat'],time_zone);

% define start date
start_date = dateshift(edf_files.start(1),'start','day');

% define end date
end_date = dateshift(edf_files.end(end),'end','day');

% import timestamps of sqEEG seizures into a table - import_EEG_seizures_v2 function, reads specific annotation labels
EEG_seizures = import_EEG_seizures_v2(subject,path_EEG_review,time_zone,{'sz_50' 'sz_75'});

%% Import and align hand-selected spikes, aligned for each channel

% read .txt file with a set of 10-60 spikes marked manually on EpiSight - careful to mark close to the
% sharpest negative or positive peak of the discharge
T = readtable([path_IED filesep subject filesep subject '_sp_rev.txt']);
ied_dt = {};

% spikes are labelled "sp_,_", first number for template number (if more than one template), second
% number for channel where spike is seen ("sp_" only if spike is seen on
% both channels).

ti = 1;
while any(contains(T.Var4,['sp' num2str(ti)]))

ied_dt(ti,:) = {datetime(T.Var2,'TimeZone',time_zone) datetime(T.Var2,'TimeZone',time_zone)};

ied_dt{ti,1} = ied_dt{ti,1}(strcmp(T.Var4,['sp' num2str(ti)]) | strcmp(T.Var4,['sp' num2str(ti) ',1']));
ied_dt{ti,2} = ied_dt{ti,2}(strcmp(T.Var4,['sp' num2str(ti)]) | strcmp(T.Var4,['sp' num2str(ti) ',2']));

ied_dt{ti,1}=sort(ied_dt{ti,1});
ied_dt{ti,2}=sort(ied_dt{ti,2});
% 
ti = ti+1;
end

%% Align spikes

% extract hdr of a random EEG file to obtain the sampling rate
hdr = ft_read_header('.edf');
srate = hdr.Fs;

% define peri-spike time (set to 1 second before and after the peak)
peri_spike_time = [1 1]; % in seconds
n_samples = sum(round(peri_spike_time*srate))+1;
max_lag = 20; % 20 samples before and after peak allowed

max_prom = {}; loc_max_prom = {}; max_neg_prom = {}; loc_max_neg_prom = {}; pk_lag = {};
ied_seg = {}; ied_seg_time = {};

load gong

% loop over templates, channels and spikes
figure
t = tiledlayout('flow');

for iti=1:2 % two iterations, one for alignment and the other to retrieve spikes
    for ti = 1:size(ied_dt,1)
        for chi = 1:2
            nexttile
            if iti==1
                title(['Temp. ' num2str(ti) ', chan. ' num2str(chi) ' - misaligned'])
            else
                title(['Temp. ' num2str(ti) ', chan. ' num2str(chi) ' - aligned'])
            end
            hold on
            for i = 1:size(ied_dt{ti,chi},1)
                disp(['Template ' num2str(ti) ', channel ' num2str(chi) ', spike ' num2str(i)])

                % find which file does the spike belong to
                A = ied_dt{ti,chi}(i) - edf_files.start;
                ied_file(i) = find(A == min(A(A>0)));

                % read files (skip reading if same file has been loaded)
                if i<2 || (i>1 && ied_file(i) ~= ied_file(i-1))

                    load([path_EEG filesep subject filesep 'mat' filesep char(edf_files.filename(ied_file(i)))]);

                % define time vector
                time = edf_files.start(ied_file(i)):seconds(1/srate):edf_files.end(ied_file(i));
                time = time(1:end-1);
                else, end

                % select only eeg/time containing spike (+- peri-spike time)
                [a idx] = min(abs(ied_dt{ti,chi}(i) - time));

                ied_seg{ti,chi}(i,:) = data(chi,idx-round(peri_spike_time(1)*srate):idx+round(peri_spike_time(end)*srate));
                ied_seg_time{ti,chi}(i,:) = time(idx-round(peri_spike_time(1)*srate):idx+round(peri_spike_time(end)*srate));

                plot(ied_seg{ti,chi}(i,:))


                if iti==1
                    [pks,locs,~,p] = findpeaks(ied_seg{ti,chi}(i,round(n_samples/2)-max_lag:round(n_samples/2)+max_lag));
                    if isempty(pks)
                        max_prom{ti,chi}(i) = 0; ind_pk = 0; loc_max_prom{ti,chi}(i) = 0;
                    else
                        [max_prom{ti,chi}(i),ind_pk] = max(p);
                        loc_max_prom{ti,chi}(i) = locs(ind_pk);
                    end

                    [pks_neg,locs_neg,~,p] = findpeaks(-ied_seg{ti,chi}(i,round(n_samples/2)-max_lag:round(n_samples/2)+max_lag));
                    if isempty(pks_neg)
                        max_neg_prom{ti,chi}(i) = 0; ind_neg_pk = 0; loc_max_neg_prom{ti,chi}(i) = 0;
                    else
                        [max_neg_prom{ti,chi}(i),ind_neg_pk] = max(p);
                        loc_max_neg_prom{ti,chi}(i) = locs_neg(ind_neg_pk);
                    end
                end


            end
            sound(y,Fs)
            % Spike alignment
            if iti==1
                [~,b] = max([max_prom{ti,chi} ; max_neg_prom{ti,chi}]);
             
                if median(b) == 1
                    pk_lag{ti,chi} = floor(n_samples/2)-max_lag+loc_max_prom{ti,chi};
                else
                    pk_lag{ti,chi} = floor(n_samples/2)-max_lag+loc_max_neg_prom{ti,chi};
                end

                for spi=1:size(ied_dt{ti,chi},1)
                    ied_dt{ti,chi}(spi) = ied_seg_time{ti,chi}(spi,pk_lag{ti,chi}(spi));
                end
            end
        end
    end
end

save([path_IED filesep subject filesep 'first_templates.mat'],'ied_seg*','srate')
saveas(gcf,[path_IED filesep subject filesep 'spike_alignment.png'])

%% Plot spikes and averages

time2plot = (-peri_spike_time(1):1/srate:peri_spike_time(2));  

for ti = 1:size(ied_seg,1)
    for chi = 1:2
        figure
        t = tiledlayout('flow');
        for i = 1:size(ied_seg{ti,chi},1)
            nexttile
            title(num2str(i))
            data2plot = ied_seg{ti,chi}(i,:);
            plot(time2plot,data2plot)
            datetick('x','HH:MM:ss')
            xticklabels({})
            set(gca,'XLim',[time2plot(1) time2plot(end)])
            set(gca,'YLim',[-max(max(max(abs(ied_seg{ti,chi})))) max(max(max(abs(ied_seg{ti,chi}))))])
            hold on
            plot([peri_spike_time(1) peri_spike_time(1)],get(gca,'YLim'),'k-','LineWidth',0.7)
        end
        title(t,['Manual spikes - template ' num2str(ti) ', channel ' channel_label{chi}]);
        t.Padding = 'compact';
        t.TileSpacing = 'compact';
        saveas(gcf,[path_IED filesep subject filesep get(get(t,'Title'),'string') '.png'])
        
        figure
        for i = 1:size(ied_seg{ti,chi},1)
            plot(time2plot,ied_seg{ti,chi}(i,:))
            hold on
        end
        plot(time2plot,mean(ied_seg{ti,chi}),'k-','LineWidth',3)
        set(gca,'XLim',[time2plot(1) time2plot(end)]);
        xlabel('Time (sec.)')
        ylabel('\muV')
        title(['Template ' num2str(ti) ', channel ' channel_label{chi} ' - ' num2str(size(ied_seg{ti,chi},1)) ' spikes.'])
        saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
    end
end

%% Select short ied segment based on visual analysis and extract first template(s)

% define onset and offset of spike by looking at the plot, or selecting the same for all spikes
for ti=1:size(ied_seg,1)
%     prompt = ['Input spike_onset (in sec.) from template ' num2str(ti) ':'];
%     spike_onset(ti) = input(prompt);
spike_onset(ti) = -0.1;

%     prompt = ['Input spike_offset (in sec.) from template ' num2str(ti) ':'];
%     spike_offset(ti) = input(prompt);
spike_offset(ti) = 0.2;

idx = dsearchn(time2plot',spike_onset(ti)):dsearchn(time2plot',spike_offset(ti));

% Prepare first template
figure
tiledlayout('flow')
for ti = 1:size(ied_seg,1)
    for chi = 1:2
        
        ied_seg_short{ti,chi} = ied_seg{ti,chi}(:,idx);
        ied_seg_time_short{ti,chi} = ied_seg_time{ti,chi}(:,idx);
        
        % identify max and min. peak to peak from first template spikes
        maximum_p2p(ti,chi) = ceil(max(max(ied_seg_short{ti,chi},[],2) - min(ied_seg_short{ti,chi},[],2)));
        minimum_p2p(ti,chi) = floor(min(max(ied_seg_short{ti,chi},[],2) - min(ied_seg_short{ti,chi},[],2)));

        % z-score
        ied_seg_short_z{ti,chi} = (ied_seg_short{ti,chi} - mean(ied_seg_short{ti,chi},2)) ./ std(ied_seg_short{ti,chi},[],2);
        
        % first template
        templates{ti,chi} = mean(ied_seg_short_z{ti,chi});

        time2plot_short = (0:size(ied_seg_time_short{ti,chi},2)-1)./srate;

        nexttile
        for i = 1:size(ied_seg_short_z{ti,chi},1)
            plot(time2plot_short,ied_seg_short_z{ti,chi}(i,:))
            hold on
        end
        plot(time2plot_short,templates{ti,chi},'k-','LineWidth',5)
        title(['Template ' num2str(ti) ', channel ' channel_label{chi}])
        ylabel('z-score'), xlabel('Time (sec.)')
        set(gca,'XLim',[time2plot_short(1) time2plot_short(end)])
%         saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
    end
end
end

save([path_IED filesep subject filesep 'first_templates.mat'],'templates','ied_seg*','srate','spike_o*','maximum_p2p','minimum_p2p')


%% Run the first template on data
% run the first template on the whole or part of dataset

minimum_r = 0.9;

  spikes_table_screen = spike_data_corr(templates,minimum_r,minimum_p2p,maximum_p2p,edf_files,[path_EEG filesep subject filesep 'mat'],subject,srate);
    % already removes spikes with overlap (second function)

%% Select minimum number 200 spikes, with highest r value possible

rs = 0.9:0.01:0.98;
n_spikes=[];
for ri = 1:length(rs)
    n_spikes(ri) = sum(spikes_table_screen{ti}.r > rs(ri))  
end
if any(n_spikes>200)
rs = rs(n_spikes>200);
else
    rs = rs(1);
end
min_r = rs(end);
spikes_table_screen{ti} = spikes_table_screen{ti}(spikes_table_screen{ti}.r > min_r,:);

%% extract final template

for ti = 1:size(templates,1)
[ied_seg_2{ti} ied_seg_time_2{ti}] = import_EEG_segments([path_EEG filesep subject filesep 'mat'],edf_files,spikes_table_screen{ti}.dt,0,spike_offset-spike_onset,srate);
end

%% Plot distribution of highly correlated spikes

% can change these values and see differences in distribution
upper_thresh = maximum_p2p; 
lower_thresh = minimum_p2p;
upper_thresh = Inf; 
lower_thresh = 0;
% min_r = 0.95;

% plot
figure
x1 = dateshift(spikes_table_screen{ti}.dt(spikes_table_screen{ti}.r > min_r & spikes_table_screen{ti}.p2p<upper_thresh & spikes_table_screen{ti}.p2p>lower_thresh),'start','day');
y1 = timeofday(spikes_table_screen{ti}.dt(spikes_table_screen{ti}.r > min_r & spikes_table_screen{ti}.p2p<upper_thresh & spikes_table_screen{ti}.p2p>lower_thresh));
plot(x1,y1,'b.')
hold on
% plot(x1,y1+1,'b.')
set(gca,'XLim',[x1(1)-hours(12) x1(end)+hours(12)])
title('Distribution of final template spikes.')
saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
savefig([path_IED filesep subject filesep get(get(gca,'Title'),'string') '.fig'])


%% Plot spikes and averages

plot_spikes(ied_seg_time_2{ti},ied_seg_2{ti},spikes_table_screen{ti}.chan,srate,subject,path_IED)
saveas(gcf,[path_IED filesep subject filesep 'final template spikes.png'])
% Adjust maximum and minimum p2p
% Check plotted spikes
% E02 - minimum p2p acceptable (20uV); maximump2p -- larger than 68.6699
% maximum_p2p = 69;

%% Final template

for ti = 1:size(templates,1)
for chi = 1:2
final_template{ti,chi} = mean(zscore(ied_seg_2{ti}(chi,:,spikes_table_screen{ti}.chan==chi),[],2),3);
end
end

save([path_IED filesep subject filesep 'final_template.mat'],'final_template','ied_seg_2', 'minimum_p2p', 'maximum_p2p','ied_seg_time_2','spikes_table_screen')

%% Validation of the final template - selecting the optimal threshold
% load([path_IED filesep subject filesep 'random_epochs.mat'])
% concat_data = concatenate_segments([path_EEG filesep 'mat'],randomepochs,edf_files,srate);

[r_final correl_m_a pval] = ROC_analysis_sqEEG(subject,path_EEG,path_IED,edf_files,time_zone,srate,EEG_seizures)

%% FINAL SPIKES 

spikes_table_final = spike_data_corr_full(templates,r_final,minimum_p2p,maximum_p2p,edf_files,[path_EEG filesep subject filesep 'mat'],subject,srate);

save([path_IED filesep subject filesep 'spike_data.mat'],'final_template','spikes_table_final','r_final', 'minimum_p2p', 'maximum_p2p')
end

%% Figure 2 - All templates
ti=1

figure
t=tiledlayout(7,2)
for subji=1:7
    
    subject=subjects_IED{1,subji}
load([path_IED filesep subject filesep 'final_template.mat'])
ied_seg_2{ti} = (ied_seg_2{ti} - mean(ied_seg_2{ti},2)) ./ std(ied_seg_2{ti},[],2);
for chi=1:2
    nexttile
    
for i = 1:sum(spikes_table_screen{ti}.chan == chi)
    plot((1:size(ied_seg_time_2{ti},2))./srate,ied_seg_2{ti}(chi,:,max(find(spikes_table_screen{ti}.chan == chi,i,'first'))),'Color',[161 197 203]./256,'LineWidth',0.5)
    hold on
end
plot((1:size(ied_seg_time_2{ti},2))./srate,mean(ied_seg_2{ti}(chi,:,find(spikes_table_screen{ti}.chan == chi)),3),'k-','LineWidth',2);
if chi==1, ylabel(subject), end
set(gca,'XLim',[1/srate size(ied_seg_time_2{ti},2)/srate]);
if subji~=7
set(gca,'XTick',[])
end
set(gca,'YTick',[])
% % set(gca,'Ylim',[-150 150])
% xlabel('Time (sec.)')
% ylabel('z-score')
if subji==1, title(channel_label{chi}), end
if subji==7, xlabel('Time (sec.)'), end
% title([subject '-' channel_label{chi}]) 
saveas(gcf,[path_IED filesep 'all_templates.png'])

end
end

t.TileSpacing = 'compact'

print('-dtiff','-r600',[path_IED filesep 'alltemplates'])

%% Figure 4 - Spike rate over time (SUB002)

subjects_IED = {'E02','E04','E06','E08','E09','SUB001','SUB002'; 'Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/Copenhagen','Europe/London','Europe/London'}
for subji = [1 3 4 5 7] % 1:length(subjects_IED)
    close all
subject = subjects_IED{1,subji}

% timezone
time_zone = subjects_IED{2,subji}

% path to diary data
path_diary = 'G:\Other computers\My MacBook Pro\Documents\Neurologia\2019\SUBER\Analysis\Diary_Data';

% path to EEG review data
path_EEG_review = 'G:\Other computers\My MacBook Pro\Documents\Neurologia\2019\SUBER\Analysis\EEG_Review';

% path to EEG data
path_EEG = 'C:\Users\Widex\Documents\SUBER\Data';
% path_EEG = '/Users/pedrofaroviana/OneDrive - King''s College London/Documents/EpiSight/EDFtestFiles/E04/EDF_E04/E04 with Implant ID'
% path_EEG = '/Users/pedrofaroviana/OneDrive - King''s College London/Documents/EpiSight/EDFtestFiles/E06'
% path to ied detection analysis
path_IED = 'G:\Other computers\My MacBook Pro\Documents\Neurologia\2019\SUBER\Analysis\Spike_detection';

% channel labels
channel_label = {'D-C','P-C'};

% read headers, onset and offset times of files - adherence_table function
edf_files = adherence_table([path_EEG filesep subject filesep 'mat'],time_zone);

% define start date
start_date = dateshift(edf_files.start(1),'start','day');

% define end date
end_date = dateshift(edf_files.end(end),'end','day');

% import times of sqEEG seizure occurrence
EEG_seizures = import_EEG_seizures_v2(subject,path_EEG_review,time_zone,{'sz_50' 'sz_75'});

load([path_IED filesep subject filesep 'spike_data.mat'])

% Spike rate

whole_time = (dateshift(start_date,'start','day') + minutes(1:minutes(dateshift(end_date,'start','day') - dateshift(start_date,'start','day'))))';
spikes_final_min = dateshift(spikes_table_final{ti}.dt,'start','minute');

if any(diff(isdst(whole_time)) == 1)
    dst_idx = find(diff(isdst(whole_time)) == 1);
    whole_time = [whole_time(1:dst_idx); NaT(60,1,'TimeZone',time_zone); whole_time(dst_idx+1:end)];
end

if any(diff(isdst(whole_time)) == -1)
    dst_idx = find(diff(isdst(whole_time)) == -1);
    whole_time = [whole_time(1:dst_idx); whole_time(dst_idx+61:end)];
end

spikes_idx = dsearchn( datenum(whole_time),datenum(spikes_final_min));
spike_rate = [(accumarray(spikes_idx,1)) ; zeros(length(whole_time)-(spikes_idx(end)),1)]; 

figure, plot(whole_time,spike_rate,'k-')
title(['|r_x_y| > ' num2str(r_final) ' - mean spikes/min: ' num2str(round(mean(spike_rate,'omitnan'),2))])
ylabel('Spikes / min')
datetick('x','dd/mm')
savefig(gcf,['Spike rate_min - final r' num2str(round(mean(spike_rate),2)) '.fig'])

spike_rate_mat = reshape(spike_rate,1440,[]);
whole_time_mat = reshape(whole_time,1440,[]);

% Plot adherence matrix
% whole_time = datetime(whole_time,'ConvertFrom','datenum','TimeZone','Europe/Lisbon');

adherence_vec = zeros(size(whole_time));

for i = 1:length(whole_time)
adherence_vec(i)= any(isbetween(whole_time(i),edf_files.start,edf_files.end));
end
adherence_mat = reshape(adherence_vec,1440,[]);
spike_rate_mat(adherence_mat == 0) = nan;
spike_rate(adherence_vec==0) = nan;

spike_rate_z = (spike_rate - mean(spike_rate,'omitnan')) / std(spike_rate,'omitnan');
% spike_rate_z = smoothdata(spike_rate_z,60,1);
spike_rate_z = smoothdata(spike_rate_z,1,"gaussian",30,"omitnan")
spike_rate_mat_z = reshape(spike_rate_z,1440,[]);

x = datenum(dateshift(EEG_seizures.onset,'start','day'));
y = datenum(EEG_seizures.onset-floor(days(EEG_seizures.onset-start_date)));


% x1 = datenum(dateshift(diary_seizures_dt,'start','day'));
% y1 = datenum(diary_seizures_dt-floor(days(diary_seizures_dt-start_date)));

spike_rate_norm = (spike_rate-min(spike_rate,[],1))./(max(spike_rate,[],1)-min(spike_rate,[],1));
spike_rate_norm_mat = reshape(spike_rate_norm,1440,[]);

figure
clf

t=tiledlayout(3,6)

nexttile([2 4])
spike_rate_norm_mat_sm = smoothdata(spike_rate_norm_mat,1,"gaussian",30,'omitnan');
imagesc(1:size(adherence_mat,2),flipud((1:1440))./1440*24,spike_rate_norm_mat_sm,'AlphaData',~isnan(spike_rate_mat))
% imagesc(whole_time_mat(1,:)',whole_time_mat(:,1),spike_rate_mat_z,'AlphaData',~isnan(spike_rate_mat))
% datetick('y','HH:MM')
% datetick('x','dd/mm')
addpath('C:\Users\Widex\Downloads')
c=makeColorMap([0 118 192]./256,[171 180 125]./256,[255 223 118]./256,50);
colormap(c)
set(gca,'YTick',[0 6 12 18 24],'YTickLabel',{'12PM','6AM','12AM','6PM','12PM'})
hold on
x=days(dateshift(EEG_seizures.onset,'start','day')-whole_time(1))+1;
y=datenum(timeofday(dateshift(EEG_seizures.onset,'start','minute'))).*24;

plot(x,y,'*','Color',[163 2 52]./265,'LineWidth',1.2,'MarkerSize',10,'DisplayName','EEG seizures')
% plot(x1,y1,'gd','LineWidth',3,'MarkerSize',10,'DisplayName','Reported seizures')
legend
set(gca,'XLim',[0.5 size(adherence_vec,1)/1440+0.5])
% title('Spike Rate (Z score) throughout study')
colorbar
saveas(gcf,[path_IED filesep subject filesep 'studyvstimeday_spike rate.png'])

nexttile([2 2])

h1 = polarhistogram(2*pi*datenum(timeofday(EEG_seizures.onset)),12,'Normalization','probability','FaceColor',[163 2 52]./256,'LineStyle','none')    
% _angle(:,idx),'FaceColor',[0 118 192]./256,'FaceAlpha',0.5,'Normalization','count','LineWidth',0.1) % edit number of bins if necessary

ax = gca;
ax.RTick = [];
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
% ax.RLim = [0 1]
% ax.RTick = [0 0.25 0.5 0.75 1]
hold on 


spike_rate_tday = mean(spike_rate_norm_mat,2,'omitnan');
spike_rate_tday = [spike_rate_tday;spike_rate_tday;spike_rate_tday]
spike_rate_tday_sm = smoothdata(spike_rate_tday,1,"gaussian",60,'omitnan')
spike_rate_tday_sm = spike_rate_tday_sm(1441:2880);
polarplot(2*pi*(1:1440)./1440,spike_rate_tday_sm,'k-','LineWidth',1)

ax.RTickLabel = [];
% set(gca,'FontSize',8)
% ax.ThetaTick = linspace(0,360-(360/),6)
% % ax.RTick = [];
% ax.FontName = 'Calibri';
% if danish(pt).cycle_lengths(idx) == 1
h1.BinEdges = 0:2*pi/12:2*pi;
ax.ThetaTickLabel = {'12am';'2am';'4am';'6am';'8am';'10am';'12pm';'2pm';'4pm';'6pm';'8pm';'10pm'};
% ax.ThetaTickLabel = {'12am';'4am';'8am';'12pm';'4pm';'8pm'};
ax.ThetaTick

nexttile([1 4])

% spike_rate_z = (spike_rate - mean(spike_rate,'omitnan')) / std(spike_rate,'omitnan');
% spike_rate_z = smoothdata(spike_rate_z,60,1);
spike_rate_norm_sm = smoothdata(spike_rate_norm,1,"gaussian",1440*7,"omitnan");

plot((1:size(adherence_vec,1))./1440,spike_rate_norm_sm,'k','LineWidth',1)

x=days(EEG_seizures.onset-whole_time(1))+1;

y = spike_rate_norm_sm(dsearchn(((1:size(adherence_vec,1))./1440)',x));
hold on, 
plot(x,y,'*','Color',[163 2 52]./265,'LineWidth',1.2,'MarkerSize',10,'DisplayName','EEG seizures')

set(gca,'XLim',[0.5 size(adherence_vec,1)/1440+0.5])
xlabel('Days of study')


nexttile([1 2])

plot((1:1440)./60,smoothdata(spike_rate_norm_mat(:,25),1,"gaussian",60,'omitnan'),'k','LineWidth',1)
set(gca,'Xlim',[0 12])
set(gca,'XTickLabel',{'12am';'2am';'4am';'6am';'8am';'10am';'12pm'})
% saveas(gcf,[path_IED filesep subject filesep 'ultradian_spike_rate.png'])



saveas(gcf,[path_IED filesep subject filesep 'Figure 4.png'])

end
