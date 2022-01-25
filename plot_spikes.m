%% Plot spikes and averages

function plot_spikes(ied_segment_time,ied_segment,chan,srate,subject,path_IED)

channel_label = {'D-C','P-C'};
% path_IED = '/Users/pedrofaroviana/Documents/Neurologia/2019/SUBER/Analysis/Spike_detection';

% raw vals and averages
figure
t=tiledlayout('flow')
for chi = 1:2
nexttile
for i = 1:sum(chan == chi)
    plot((1:size(ied_segment_time,2))./srate,ied_segment(chi,:,max(find(chan == chi,i,'first'))))
    hold on
end
plot((1:size(ied_segment_time,2))./srate,mean(ied_segment(chi,:,find(chan == chi)),3),'k-','LineWidth',3);

set(gca,'XLim',[1/srate size(ied_segment_time,2)/srate]);
% set(gca,'Ylim',[-350 350])
xlabel('Time (sec.)')
ylabel('\muV')
title(['Channel ' channel_label{chi} ' - ' num2str(sum(chan == chi)) ' spikes.']) 
end

% z-scoring

ied_segment = (ied_segment - mean(ied_segment,2)) ./ std(ied_segment,[],2);

for chi = 1:2
nexttile
for i = 1:sum(chan == chi)
    plot((1:size(ied_segment_time,2))./srate,ied_segment(chi,:,max(find(chan == chi,i,'first'))))
    hold on
end
plot((1:size(ied_segment_time,2))./srate,mean(ied_segment(chi,:,find(chan == chi)),3),'k-','LineWidth',3);

set(gca,'XLim',[1/srate size(ied_segment_time,2)/srate]);
% set(gca,'Ylim',[-150 150])
xlabel('Time (sec.)')
ylabel('z-score')
title(['Channel ' channel_label{chi} ' - ' num2str(sum(chan == chi)) ' spikes_zscore.']) 
% saveas(gcf,[path_IED filesep subject filesep get(get(gca,'Title'),'string') '.png'])
end


