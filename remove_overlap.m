%% remove overlapping spikes (i.e. those with more than half the samples)
% inputs: spikes table; template; sampling rate

function [spikes_table] = remove_overlap(spikes_table,template,srate)
disp('Removing overlapping spikes...')
% tic
spikes_table = sortrows(spikes_table,2,'descend');
overlap_length = size(template,2)/srate/2;

overlap_length_d = overlap_length / (60*60*24);

% work with datenum, it's faster
spikes_dn = datenum(spikes_table.dt);

% define temp variable and sort it with time
spikes_dn_temp = spikes_dn;
spikes_dn_temp = sort(spikes_dn_temp);

i=1;
% while still any spikes with minimum diff, continue
while any(abs(diff(spikes_dn_temp)) < overlap_length_d)
    % reordered temp variable
    spikes_dn_temp = spikes_dn;
    idx = find(abs(spikes_dn(i) - spikes_dn) < overlap_length_d);
    % remove for both variables
    spikes_dn(idx(idx>i)) = [];
    spikes_dn_temp(idx(idx>i)) = [];
    i=i+1;
    % sort it to time and rerun the while loop
    spikes_dn_temp = sort(spikes_dn_temp);
end

spikes_table = spikes_table(dsearchn(datenum(spikes_table.dt),spikes_dn),:);

spikes_table = sortrows(spikes_table,1);
% toc
end
