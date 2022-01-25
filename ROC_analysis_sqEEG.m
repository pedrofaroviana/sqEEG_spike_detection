%% Correlation + ROC analysis of manually-reviewed data
% concat_data - chan X time X epoch_number
% manual_spikes - spikes identified manually (in datetime format)
% randomepochs - datetime of each random epoch
function [rho_final correl_m_a pval] = ROC_analysis_sqEEG(subject,path_EEG,path_IED,edf_files,time_zone,srate,EEG_seizures)
% just one epoch
randomepochs = readtable([path_IED filesep subject filesep 'random_epochs.xls']);
% randomepochs = randomepochs(1,:)
randomepochs.Var1 = datetime (randomepochs.random_beg_epochs,'TimeZone',time_zone);
randomepochs.Var2 = randomepochs.Var1 + hours(1);
randomepochs = removevars(randomepochs,{'random_beg_epochs','random_end_epochs'});
save('random_epochs.mat',"randomepochs")
% load('random_epochs.mat')

% manual spikes
manual_spikes = readtable([path_IED filesep subject filesep subject '_sp_rev.txt']);
% manual_spikes = readtable('sp_rev_15092020.txt');
manual_spikes = manual_spikes (strcmp(manual_spikes.Var4,'sp rev'),:);
manual_spikes = datetime(manual_spikes.Var2,'TimeZone',time_zone);

load([path_IED filesep subject filesep 'final_template.mat'])

concat_data = concatenate_segments([path_EEG filesep subject filesep 'mat'],randomepochs,edf_files,srate);

%% Windowed correlation
% set minimum rho
disp('correlation with validated dataset...')
minimum_rho = 0.5;

spikes_table_val = table(NaT(0,'TimeZone','Europe/Lisbon'),[],[],[],'VariableNames',{'dt','rho','chan','p2p'});
n_templates = size(final_template,1);
for si = 1:size(concat_data,3)
    tic
    disp(['Segment n. ' num2str(si) '/' num2str(size(concat_data,3))])


    ied_data = squeeze(concat_data(:,:,si));

        for ti=1:n_templates
        rho = zeros(length(ied_data)-length(final_template{ti,1})+1,2);
        peak_2_peak = zeros(length(ied_data)-length(final_template{ti,1})+1,2);
        
        % loop over channels
        for chi = 1:2
            % matrix with sliding window segments
            
            data_corr = buffer(ied_data(chi,:)',length(final_template{ti,1}),length(final_template{ti,1})-1,'nodelay');
            % correlation coef.
            rho(:,chi) = corr(data_corr,(final_template{ti,chi})');
            
            % amplitude
            peak_2_peak(:,chi) = max(data_corr) - min(data_corr);
           
             % indices for segments with > mininum rho and amp. thresholds
            high_corr_idx = find(rho(:,chi)>minimum_rho & peak_2_peak(:,chi) > minimum_p2p(ti,chi) & peak_2_peak(:,chi) < maximum_p2p(ti,chi));
            spikes_table_t_chan{ti,chi} = table((randomepochs.Var1(si) + seconds((high_corr_idx./srate))), rho(high_corr_idx,chi), (zeros(length(high_corr_idx),1) + chi), peak_2_peak(high_corr_idx,chi));

        end

spikes_table_t{ti,si} = cat(1,spikes_table_t_chan{ti,:});
spikes_table_t{ti,si}.Properties.VariableNames = {'dt','rho','chan','p2p'};
spikes_table_t{ti,si} = remove_overlap(spikes_table_t{ti,si},final_template{ti,1},srate);


        end
end

spikes_table_val = cat(1,spikes_table_t{ti,:});

% % matrify computation
% data_corr = zeros(2,length(final_template),length(ied_data)-length(final_template)+1);
% 
% spikes_table_t_chan = {[] []};
% rho = zeros(length(ied_data)-length(final_template)+1,2);
% peak_2_peak = zeros(length(ied_data)-length(final_template)+1,2);
% for chi = 1:2
%     squeeze_data = ied_data(chi,:);
% data_corr(chi,:,:) = (squeeze_data((bsxfun(@plus,(1:length(final_template)).',0:1:numel(squeeze_data)-length(final_template)))'))';
% rho(:,chi) = corr(squeeze(data_corr(chi,:,:)),(final_template(chi,:))');
% peak_2_peak(:,chi) = squeeze(max(data_corr(chi,:,:)) - min(data_corr(chi,:,:)));
% high_corr_idx = find(rho(:,chi)>minimum_rho);
% spikes_table_t_chan{chi} = table((randomepochs.Var1(si) + seconds((high_corr_idx./srate))) , rho(rho(:,chi)>minimum_rho,chi), (zeros(length(high_corr_idx),1) + chi),peak_2_peak(rho(:,chi)>minimum_rho,chi));
% end
% 
% spikes_table_temp = [spikes_table_t_chan{1};spikes_table_t_chan{2}];
% spikes_table_temp.Properties.VariableNames = {'dt','rho','chan','p2p'};
% spikes_table_temp = remove_overlap(spikes_table_temp,final_template,srate);
% 
% spikes_table_val(end+1:end+size(spikes_table_temp,1),:) = spikes_table_temp;
% 
% clear rho* data_corr spikes_table_temp
% toc
% end

% save([path_IED filesep subject filesep 'spikes_table_val.mat'],'spikes_table_val')

%% Overall ROC curve
disp('Sensitivity/Specificity analysis...')
% Remove spikes during ictal periods and time periods without data
ictal_time_2_discard = 0;
non_rec_time = 0;
for i = 1:length(EEG_seizures.onset)
    if any(isbetween(EEG_seizures.onset(i), randomepochs.Var1, randomepochs.Var2))
        ictal_time_2_discard = ictal_time_2_discard + seconds(EEG_seizures.offset(i)-EEG_seizures.onset(i));
    end
spikes_table_val(isbetween(spikes_table_val.dt,EEG_seizures.onset(i),EEG_seizures.offset(i)),:) = [];
manual_spikes(isbetween(manual_spikes,EEG_seizures.onset(i),EEG_seizures.offset(i))) = [];
end

% remove time periods without data
time_2_discard = ictal_time_2_discard + ((sum(sum(sum(isnan(concat_data)))) / 2) / srate);

% Define peak to peak of template
% peak_2_peak = squeeze(max(ied_seg_2{ti},[],2) - min(ied_seg_2{ti},[],2));
% peak_2_peak = peak_2_peak(:);

% Define threshold values for p2p and rho minimum values
% threshold = [0 (0.95*prctile(peak_2_peak,0)) prctile(peak_2_peak,0:5)];
rho = 0.5:0.01:1;

% Define positives and negatives
num_manual_spikes_total = length(manual_spikes);
num_negatives = round(((size(concat_data,2)/srate)*size(concat_data,3) - time_2_discard) /0.3 - length(manual_spikes)); 

% Initialize variables
matched_spikes = cell(1,length(rho));
num_autom_spikes_total = zeros(1,length(rho));
tpositives_total = zeros(1,length(rho));
fpositives_total = zeros(1,length(rho));

% matched_spikes = cell(length(threshold),length(rho));
% num_autom_spikes_total = zeros(length(threshold),length(rho));
% tpositives_total = zeros(length(threshold),length(rho));
% fpositives_total = zeros(length(threshold),length(rho));


% for thri = 1:length(threshold)
thri=1;
for rhi = 1:length(rho)

manual_spikes_t2 = manual_spikes;

autom_spikes_total = spikes_table_val.dt(spikes_table_val.rho > rho(rhi));

correspondence = zeros(size(autom_spikes_total));
% matched_spikes = [];

for i = 1:length(autom_spikes_total)
if min(abs(seconds(autom_spikes_total(i) - manual_spikes_t2))) < 0.3 % any sample overlaps in time
[a b] = min(abs(seconds(autom_spikes_total(i) - manual_spikes_t2)));
correspondence(i) = 1;
matched_spikes{thri,rhi} = [matched_spikes{thri,rhi} ;manual_spikes_t2(b)];
manual_spikes_t2(b) = []; 
else, end
end

num_autom_spikes_total(thri,rhi) = length(autom_spikes_total);

tpositives_total(thri,rhi) = length(matched_spikes{thri,rhi});
fpositives_total(thri,rhi) = sum(correspondence == 0);

end

% end

tpositiverate_total = tpositives_total ./ num_manual_spikes_total;

fpositiverate_total = fpositives_total ./ num_negatives;

ratio_fpositive_spikes = fpositives_total ./ num_autom_spikes_total;
ratio_fpositive_spikes(isnan(ratio_fpositive_spikes)) = 0;

fpositiverate_hour_total = fpositives_total ./ size(concat_data,3);

percent_rate_total = num_autom_spikes_total ./ length(manual_spikes);

%% Plot ROC Curves
disp('plotting ROC curves...')
% Ratio of false positive spikes
figure
t=tiledlayout('flow');
nexttile
plot_cmap = [zeros(size(fpositiverate_total,1),1) zeros(size(fpositiverate_total,1),1) (linspace(0,1,size(fpositiverate_total,1)))'];
for i = 1:size(fpositiverate_total,1)
    plot(ratio_fpositive_spikes(i,:),tpositiverate_total(i,:),'Color',plot_cmap(i,:),'LineWidth',2)
    hold on
end
% legend([string(threshold),"random"])
hold on, plot([0 1],[0 1],'k--')
set(gca,'XLim',[0 1]), set(gca,'YLim',[0 1])
axis square
title(['ROC curves - ratio of positive spikes to total detected spikes'])
xlabel('Ratio of false positive spikes'), ylabel('True Positive Rate')
% saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])

% "Real" false positive rate
nexttile
plot_cmap = [zeros(size(fpositiverate_total,1),1) zeros(size(fpositiverate_total,1),1) (linspace(0,1,size(fpositiverate_total,1)))'];
for i = 1:size(fpositiverate_total,1)
    plot(fpositiverate_total(i,:),tpositiverate_total(i,:),'Color',plot_cmap(i,:),'LineWidth',2)
    hold on
end
% legend([string(threshold),"random"])
hold on, plot([0 1],[0 1],'k--')
set(gca,'XLim',[0 1]), set(gca,'YLim',[0 1])
axis square
title(['ROC curves - real false positive rate'])
xlabel('False Positive Rate'), ylabel('True Positive Rate')
% saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])

% False alarms/hour
nexttile
plot_cmap = [zeros(size(fpositiverate_hour_total,1),1) zeros(size(fpositiverate_total,1),1) (linspace(0,1,size(fpositiverate_total,1)))'];
for i = 1:size(fpositiverate_hour_total,1)
    plot(fpositiverate_hour_total(i,:),tpositiverate_total(i,:),'Color',plot_cmap(i,:),'LineWidth',2)
    hold on
end
% legend([string(threshold),"random"])
hold on, plot([0 3600/0.15],[0 1],'k--')
set(gca,'XLim',[0 3600/0.15]), set(gca,'YLim',[0 1])
axis square
title(['ROC curves - FAR per hour'])
xlabel('False positive spikes / hour'), ylabel('True Positive Rate')

% saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
saveas(gcf,[path_IED filesep subject filesep subject '_ROC_curves.png'])

% Youden Index
youdenindex = tpositiverate_total + (1-ratio_fpositive_spikes) - 1;
[youden_max, youden_max_idx] = max(youdenindex,[],2);
rho_final = rho(youden_max_idx)';

% Closest to [0 1] Criteria
youdenindex2 = sqrt(ratio_fpositive_spikes.^2 + (1-tpositiverate_total).^2);
[youden_max2, youden_max_idx2] = max(youdenindex2,[],2);
rho_final2 = rho(youden_max_idx2)';

% Concordance Probability Method - selected
youdenindex = tpositiverate_total .* (1-ratio_fpositive_spikes);
[youden_max, youden_max_idx] = max(youdenindex,[],2);
rho_final = rho(youden_max_idx)';

youdenindex_real = tpositiverate_total + (1-fpositiverate_total) - 1;
[youdenindex_real_max, youden_index_real_idx] = max(youdenindex_real,[],2);


for i = 1:length(youden_max_idx)
fpositive_youden(i,:) = ratio_fpositive_spikes(i,youden_max_idx(i));
tpositive_youden(i,:) = tpositiverate_total(i,youden_max_idx(i));
fpositive_h_youden(i,:) = fpositiverate_hour_total(i,youden_max_idx(i));
end



[mval, max_idx] = max(youden_max);


% p2pfinal = max_values.p2p_threshold(max_idx);

%% Plot selected spikes
disp('Plotting data segments...')
for ri = 1:size(randomepochs,1)
youdenspikes = spikes_table_val(spikes_table_val.rho > rho_final,:);

youdenspikes = youdenspikes(isbetween(youdenspikes.dt,randomepochs.Var1(ri),randomepochs.Var2(ri)),:);
for chi = 1:2
chan_sp{chi} = youdenspikes.dt(youdenspikes.chan == chi);
end

manual_spikes_t = manual_spikes(isbetween(manual_spikes,randomepochs.Var1(ri),randomepochs.Var2(ri)));


time = datenum(randomepochs.Var1(ri) + seconds((1:length(concat_data))./srate));

figure, 
t = tiledlayout(3,1);
for chi = 1:2
    ax(chi) = nexttile;
    plot(time,squeeze(concat_data(chi,:,ri)))
    hold on
%     plot(datenum(chan_sp{chi}),concat_data(chi,dsearchn(time',datenum(chan_sp{chi})),ri),'r*','MarkerSize',10)
    plot(datenum(chan_sp{chi}),zeros(1,length(chan_sp{chi})),'r*','MarkerSize',10)
    set(gca,'YLim',[-100 100])
%     plot(datenum(manual_spikes_t),concat_data(chi,dsearchn(time',datenum(manual_spikes_t)),ri),'gd','MarkerSize',10)
    plot(datenum(manual_spikes_t),zeros(1,length(manual_spikes_t)),'gd','MarkerSize',10)
    datetick('x','HH:MM')
end
ax(3) = nexttile;
plot(datenum(youdenspikes.dt),zeros(1,length(youdenspikes.dt))-50,'r*','MarkerSize',10,'LineWidth',2)
hold on
plot(datenum(manual_spikes_t),zeros(1,length(manual_spikes_t))+50,'gd','MarkerSize',10,'LineWidth',2)
set(gca,'YLim',[-1 1]), set(gca,'YTickLabel',[]), set(gca,'YTick',[])
    datetick('x','HH:MM')
linkaxes([ax(:)])
set(gca,'XLim',[datenum(randomepochs.Var1(ri)) datenum(randomepochs.Var2(ri))])
title(t,['Random segment ' num2str(ri) ])
saveas(gcf,[path_IED filesep subject filesep get(get(t,'Title'),'string') '.png'])
end

%% Stats for each epoch (with my final rho)
disp('Stats for each epoch...')
autom_spikes = [];

num_autom_spikes = zeros(size(randomepochs,1),1);
tpositives = zeros(size(randomepochs,1),1);
num_manual_spikes = zeros(size(randomepochs,1),1);
fpositives = zeros(size(randomepochs,1),1);

for ri = 1:size(randomepochs,1)    
manual_spikes_t = manual_spikes(isbetween(manual_spikes,randomepochs.Var1(ri),randomepochs.Var2(ri)));
manual_spikes_t2 = manual_spikes_t;

autom_spikes = spikes_table_val(spikes_table_val.rho > rho_final,:);
autom_spikes = autom_spikes.dt(isbetween(autom_spikes.dt,randomepochs.Var1(ri),randomepochs.Var2(ri)),:);
correspondence = zeros(size(autom_spikes));
matched_spikes = [];
for i = 1:length(autom_spikes)
if min(abs(seconds(autom_spikes(i) - manual_spikes_t2))) < 0.3 % any sample overlaps in time
[a b] = min(abs(seconds(autom_spikes(i) - manual_spikes_t2)));
correspondence(i) = 1;
matched_spikes = [matched_spikes;manual_spikes_t2(b)];
manual_spikes_t2(b) = []; 
else, end
end

num_autom_spikes(ri,:) = length(autom_spikes);
tpositives(ri,:) = length(matched_spikes);
num_manual_spikes(ri,:) = length(manual_spikes_t);
fpositives(ri,:) = sum(correspondence == 0);
end

fpositiverate = fpositives ./ num_autom_spikes;
tpositiverate = tpositives ./ num_manual_spikes;
fpositiverate_hour = fpositives;
percent_rate = num_autom_spikes ./ num_manual_spikes;

epochs_table = table((1:size(randomepochs,1))',   randomepochs.Var1,timeofday(randomepochs.Var1),num_autom_spikes,num_manual_spikes,percent_rate,tpositiverate,fpositiverate,fpositiverate_hour,'VariableNames',...
    {'epoch_n','start_date','time_of_day','num_autom_spikes','num_manual_spikes','percent_rate','tpositiverate','fpositiverate','fpositiverate_hour'})
epochs_table = sortrows(epochs_table,3)
writetable(epochs_table,[path_IED filesep subject filesep 'IED detection validation stats.xlsx'])

%% Figure 3
figure
tiledlayout(2,2)
nexttile([1 1])
scatter(num_manual_spikes,num_autom_spikes,10,'k','filled','o')
hold on
plot([0 500],[0 500],'k--')
set(gca,'Xlim',[0 500],'YLim',[0 500])
% scatter(hours(epochs_table.time_of_day),epochs_table.fpositiverate_hour,'k.')
xlabel ('n. manual spikes')
ylabel('n. automatic spikes')

nexttile([1 1])
scatter(hours(epochs_table.time_of_day),epochs_table.percent_rate,10,'k','filled','o')

set(gca,'Xlim',[0 23],'XTick',[0 6 12 18])
set(gca,'XTickLabel',{'12pm','6am','12am','6pm'})
xlabel ('time of day')
ylabel('percent rate')
nexttile([1 2])

disp('Correlation between automatic and manual spike rate...')
autom_spike_min = datenum((dateshift(spikes_table_val.dt(spikes_table_val.rho > rho_final),'start','minute')));
manual_spikes_min = datenum((dateshift(manual_spikes,'start','minute')));

val_time = [];
for i = 1:length(randomepochs.Var1)
    val_time = [val_time datenum(dateshift(randomepochs.Var1(i),'start','minute') + minutes(0:59))];
end
val_time = val_time';
autom_spikes_idx = dsearchn( val_time,autom_spike_min);
manual_spikes_idx = dsearchn( val_time,manual_spikes_min);

autom_spike_rate = accumarray(autom_spikes_idx,1);
autom_spike_rate(end:length(val_time)) = 0;
manual_spike_rate = accumarray(manual_spikes_idx,1); 
manual_spike_rate(end:length(val_time)) = 0;

[correl_m_a pval] = corr(autom_spike_rate,manual_spike_rate);

plot(autom_spike_rate,'Color',[0 118 192]./256), hold on, plot(manual_spike_rate,'Color',[206 128 128]./256)
ylabel('Spike rate (Sp / min.)'), xlabel('Time (min.)')
legend({'Autom. detection','Visual identification'})
% text(size(randomepochs,1)*60,1,['Rho = ' num2str(correl_m_a) ', p = ' num2str(pval) '   '],'HorizontalAlignment','right','FontSize',12)
% set(gca,'XLim',[0 780])
title('Correlation between manual and automated detection')
saveas(gcf,[path_IED filesep 'Figure 3.png'])
print('-dtiff','-r600',[path_IED filesep 'correl_manual_autom'])




%% Correlation between automated spike rate and manual spike rate
disp('Correlation between automatic and manual spike rate...')
autom_spike_min = datenum((dateshift(spikes_table_val.dt(spikes_table_val.rho > rho_final),'start','minute')));
manual_spikes_min = datenum((dateshift(manual_spikes,'start','minute')));

val_time = [];
for i = 1:length(randomepochs.Var1)
    val_time = [val_time datenum(dateshift(randomepochs.Var1(i),'start','minute') + minutes(0:59))];
end
val_time = val_time';
autom_spikes_idx = dsearchn( val_time,autom_spike_min);
manual_spikes_idx = dsearchn( val_time,manual_spikes_min);

autom_spike_rate = accumarray(autom_spikes_idx,1);
autom_spike_rate(end:length(val_time)) = 0;
manual_spike_rate = accumarray(manual_spikes_idx,1); 
manual_spike_rate(end:length(val_time)) = 0;

[correl_m_a pval] = corr(autom_spike_rate,manual_spike_rate);

figure, plot(autom_spike_rate), hold on, plot(manual_spike_rate)
ylabel('Spike rate (Sp / min.)'), xlabel('Time (min.)')
% plot([(0:60:size(randomepochs,1)*60)' (0:60:size(randomepochs,1)*60)'],get(gca,'YLim'),'k-')
legend({'Autom. detection','Visual identification'})
text(size(randomepochs,1)*60,1,['Rho = ' num2str(correl_m_a) ', p = ' num2str(pval) '   '],'HorizontalAlignment','right','FontSize',12)
% set(gca,'XLim',[0 780])
title('Correlation between manual and automated detection')
saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])

max_values = table(youden_max,rho(youden_max_idx)',fpositive_youden,fpositive_h_youden,  tpositive_youden,correl_m_a,pval,'VariableNames',{'Youden_idx_max','Rho','FPR','FAR/h','TPR','Correlation','pval'})

writetable(max_values,[path_IED filesep subject filesep subject '_val_stats.xlsx'])

end

% %% ROC curve for each segment
% 
% for i = 1:size(randomepochs,1)
%     figure, plot(fpositiverate(i,:),tpositiverate(i,:))
%     xlabel('False Positive Rate'), ylabel('True Positive Rate')
%     hold on, plot([0 1],[0 1],'k-')
%     set(gca,'XLim',[0 1]), set(gca,'YLim',[0 1])
%     title(['Segment ' num2str(i) ', ' char(randomepochs.Var1(i)) ' - ' num2str(num_manual_spikes(i)) ' manual spikes'])
%     saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
% end
% 
% %% NOTES --------------------------
% %% Run final template on manually reviewed data
% 
% T = readtable([path_EEG_review filesep subject filesep subject '_annotations.xlsx']);
% manual_spikes = T.Var2(contains(T.Var4,'sp rev'));
% manual_spikes = datetime(manual_spikes,'TimeZone','Europe/Lisbon');
% clear T
% 
% randomepochs = readtable([path_EEG_review filesep subject filesep 'random_dates_27082020.txt']);
% randomepochs = randomepochs(1:2:end,:);
% randomepochs.Var2 = randomepochs.Var1 + hours(1);
% randomepochs.Var1 = datetime(randomepochs.Var1,'TimeZone','Europe/Lisbon');
% randomepochs.Var2 = datetime(randomepochs.Var2,'TimeZone','Europe/Lisbon');
% 
% concat_data = concatenate_segments(randomepochs,edf_files,srate);
% 
% save([path_IED filesep subject filesep 'concat_data.mat'],'concat_data','randomepochs')
% 
% %% Comparison with manual selection
% 
% % put spikes from both channels together
% automated_spikes = table([spikes_tabnew{1}.Var1 ; spikes_tabnew{2}.Var1],[spikes_tabnew{1}.Var2 ; spikes_tabnew{2}.Var2],[ones(size(spikes_tabnew{1},1),1); ones(size(spikes_tabnew{2},1),1)+1]);
% 
% automated_spikes = sortrows(automated_spikes,2);
% 
% % remove spike detections separated by < 5 samples
% 
% automated_spikes = addvars(automated_spikes,[10; seconds(diff(automated_spikes.Var1))]);
% automated_spikes = addvars(automated_spikes,automated_spikes.Var4<0.3); % 300 miliseconds
% transitions = diff(automated_spikes.Var5);
% if transitions(end) == 1, transitions(end) = 0; else, end
% indexes = [(find([0;transitions == 1]) - 1) (find([0;transitions == -1])-1)];
% spikes_tabnew = automated_spikes;
% 
% spikes_indexes = {};
% for i =  1:size(indexes,1)
% spikes_indexes{i} = indexes(i,1):indexes(i,2);
% end
% 
% spikes_tabnew(cell2mat(spikes_indexes),:) = [];
% 
% for i =  1:size(indexes,1)
%     [max_rho(i,1) max_rho(i,2)] = max(automated_spikes.Var2(spikes_indexes{i}));
%     spikes_tabnew = [spikes_tabnew;automated_spikes(spikes_indexes{i}(max_rho(i,2)),:)];
% end
% spikes_tabnew = sortrows(spikes_tabnew,1);
% 
% %% Read manual identified spikes
% 
% % Manual spike identification
% spikes_tab_manual = readtable('sp_rev_15092020.txt');
% spikes_tab_manual.Var2 = datetime(spikes_tab_manual.Var2,'TimeZone','Europe/Lisbon')
% 
% % Random epoch identification
% randomepochs = readtable('random_dates_27082020.txt')
% randomepochs = randomepochs(1:2:end,:);
% randomepochs.Var2 = randomepochs.Var1 + hours(1);
% randomepochs.Var1 = datetime(randomepochs.Var1,'TimeZone','Europe/Lisbon')
% randomepochs.Var2 = datetime(randomepochs.Var2,'TimeZone','Europe/Lisbon')