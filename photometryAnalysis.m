close all; clear all; clc;
%% Import data
addpath('functions/')
pathtofile = 'D:\Vale\Lab Nabavi\Doric FP\Script\Input data\'; %folder name

% filename = 'Photometry 1 0L0R terminals black loom_trial 526.doric';
filename = '2021-12-14 16-36-53_GCaMP8m 0L0R cage 2 200um recall_2nd recall.doric'; 
% filename = 'Photometry 1 1L1R_1_trial 527.doric';
% filename = 'whiter cap 1L0R terminals black loom_trial 523.doric';
% filename = 'PHotometry 2 1L1R no freezig 1st loom_did not see loom at min 13_ no reaction to last loom_trial 532.doric'; %% Very noisy

%Specify data locations; NB h5disp([pathtofile, filename]) to get structure information
ref_loc = '/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)'; %location of isosbestic signal in structure
sig_loc = '/Traces/Console/AIn-1 - Dem (AOut-2)/AIn-1 - Dem (AOut-2)'; %location of signal data in structure
ttl_loc = '/Traces/Console/DI--O-1/DI--O-1'; %TTL location
time_loc = '/Traces/Console/Time(s)/Console_time(s)'; %Time location

%Import data
ref_data = h5read([pathtofile,filename],ref_loc); %import isosbestic data
sig_data = h5read([pathtofile,filename],sig_loc); %import signal data
ttl_data = h5read([pathtofile,filename],ttl_loc); %import ttl data
time_data = h5read([pathtofile,filename],time_loc); %import signal data

%% Set Parameters
bl_gap_size = 2; %time to include before ttl pulse in plots
bl_length = 2; %number of seconds to include for baseline calculations
af_gap_size =  3; %time to include in plots after events

T = diff(time_data(1:2)); %get period
sr = 1/T; %sample rate; kinda weird number?

% ndsampfac = 1; %downsample factor
ndsampfac = 10; %downsample factor
nsr = round(sr/ndsampfac); %round downsampled samplerate
nT = 1/nsr; %new period

%% Miscellaneous
%Cosmetics
lcmap = lines(7); %line colors
cmap = createcmap(256, 'magma'); %colormap for heatmap, 'inferno', 'magma', 'plasma' and 'viridis' as inputs to change, otherwise you can look up on matlab the generic ones and call colormap(cmapname)

%Remove missing data
ttl_data(isnan(ref_data)) = [];
time_data(isnan(ref_data)) = [];
sig_data(isnan(ref_data)) = [];
ref_data(isnan(ref_data)) = [];

if sum(isnan(sig_data)) ~=0
    ttl_data(isnan(sig_data)) = [];
    time_data(isnan(sig_data)) = [];
    ref_data(isnan(sig_data)) = [];
    sig_data(isnan(sig_data)) = [];
end

%% Inspect raw data
figure(1)
subplot(2,1,1)
plot(time_data, ref_data, 'k')
xlabel('Time (s)')
ylabel('Amplitude (arb. units)')
title('Raw Isosbestic signal')
axis tight

subplot(2,1,2)
plot(time_data, sig_data, 'k')
ylabel('Amplitude (arb. units)')
hold on
yyaxis right
p1 = plot(time_data, ttl_data);
xlabel('Time (s)')
title('Raw GCaMP signal')
legend(p1,'TTL locations','location','best')
set(gca,'ytick',[]);
set(gca,'ycolor',[0 0 0])

%% Downsample signal, then smooth using LOWESS(Locally Weighted Scatterplot Smoothing)
sig_ds = arrayfun(@(i)  mean(sig_data(i:i+ndsampfac-1)), 1:ndsampfac:length(sig_data)-ndsampfac+1); %Down sample signal by local averaging
ref_ds = arrayfun(@(i)  mean(ref_data(i:i+ndsampfac-1)), 1:ndsampfac:length(ref_data)-ndsampfac+1); %Downsample isosbestic by local averaging
time_ds = time_data(1:ndsampfac:end-1);
ttl_ds = ttl_data(1:ndsampfac:end-1);

sig_ds2 = smoothdata(sig_ds, 'lowess','SmoothingFactor',0.05); %Lowess smoothing of signal, mb not needed
ref_ds2 = smoothdata(ref_ds, 'lowess','SmoothingFactor',0.05); %Lowess smoothing of signal

%% Detrend and get dFF
bls = polyfit(ref_ds, sig_ds, 1); %% Fit first order polynomial, least squares method
scaledctrl = bls(1).*ref_ds + bls(2); %% Get trend in reference
dFF = 100.*(sig_ds-scaledctrl)./scaledctrl; %% Detrend signal and normalize

figure(2)
plot(time_ds(1:size(dFF,2)), dFF)
xlabel('Time (s)')
ylabel('\Delta F/F (%)')
hold on
yyaxis right
plot(time_ds, ttl_ds)
title('Detrended and normalized GCaMP signal')


%% Extract ttl events
[~, event_lcs1] = findpeaks(diff(ttl_ds)); %Extract ttl location data
[~, event_lcs2] = findpeaks(-diff(ttl_ds)); %Extract ttl location data
event_lcs = [event_lcs1, event_lcs2] + 1; %Put in matrix
event_lgth = min(event_lcs(:, 2) - event_lcs(:, 1));
event_matrix = NaN(size(event_lcs, 1), event_lgth+bl_gap_size*nsr+af_gap_size*nsr); %Preallocate event matrix, take smallest of TTL lengths
time_event_matrix = (0:nT:size(event_matrix,2)/nsr)-bl_gap_size; %Preallocate timevector

for ival = 1:size(event_lcs, 1)
    f_ind = event_lcs(ival, 1) - bl_gap_size*nsr; %get first index
    s_ind = f_ind + size(event_matrix, 2); %get second index
    event_matrix(ival, 1:(s_ind-f_ind)) = dFF(f_ind:s_ind-1); %Extract event data +/- gap_size seconds into one matrix
end

%% "Raw" dFF plot
tvec = time_event_matrix(1:size(event_matrix,2)); %time vec
dff_std = std(event_matrix,[],1,'omitnan');
errvec = std(event_matrix,[],1,'omitnan')./sqrt(size(event_matrix, 1)); %Standard error of the mean
dff_sem = errvec;
dff_mvec = mean(event_matrix,'omitnan'); %Mean

uerr = dff_mvec+errvec; %upper errorbar
lerr = dff_mvec-errvec; %lower errorbar

mvec_der = diff(dff_mvec); %Derivative of mvec, for finding peaks
mvec_der = smoothdata(mvec_der, 'lowess','SmoothingFactor',0.5);
mval = max(mvec_der); 
[~, derlcs] = findpeaks(mvec_der,'minpeakheight',mval*0.48,'minpeakdistance',2*nsr); %Find where slopes are steepest

ptmat = NaN(length(derlcs),2); %preallocate matrix for points
trnpts = mvec_der == abs(mvec_der); %logical vector to stepwise find indices
    
for ival = 1:length(derlcs)
    fwdextr = trnpts(derlcs(ival):derlcs(ival)+nsr); %Find forward extrema; MAYBE CHANGE HERE TO JUST GET PEAK
    fwdextr = min(find(fwdextr == 0)); %Find forward extrema
    
    revextr = fliplr(trnpts(derlcs(ival)-nsr:derlcs(ival))); %Find reverse extrema
    revextr = min(find(revextr == 0)); %Find reverse extrema
    
    ptmat(ival, 1) = derlcs(ival)-revextr; %Put into point matrix
    ptmat(ival, 2) = derlcs(ival)+fwdextr;
    
end

ptdiff = dff_mvec(ptmat(:,2))-dff_mvec(ptmat(:,1));
ptmat(ptdiff<1,:) = []; %remove false pts

figure(3)
h1 = subplot(211);
orgSize1 = get(gca, 'Position');
p2 = plot(tvec, event_matrix','Color',[0.7 0.7 0.7]);
hold on
p1 = patch([tvec, fliplr(tvec)], [lerr, fliplr(uerr)], lcmap(1, :), 'edgecolor','none');
yline(0,':')
p3 = xline(0, '--');
xline(event_lgth*nT, '--')
p4 = plot(tvec, dff_mvec, 'color', 'k','linewidth', 3);
p5 = plot(tvec(ptmat(:)), dff_mvec(ptmat(:)), 'Color',lcmap(3, :),'Linestyle','none','Marker','x', 'markersize', 25);
xlim([tvec(1),  tvec(end)])
xlabel('Time (s)')
ylabel('\Delta F/F (%)')
legend([p3, p2(1), p4, p1, p5], {'TTL on/offset','Trials', 'Mean response', 'SEM', 'Peaks/valleys'},'location','northwest')

h2 = subplot(212);
imagesc(tvec, 1, event_matrix)
xlabel('Time (s)')
set(gca,'YTick', 1:size(event_matrix,2))
ylabel('Trial #')
colormap(cmap)
orgSize2 = get(gca, 'Position');
cb = colorbar;
set(h1, 'Position', orgSize1)
set(h2,'Position', orgSize2)
cb.Label.String = '\DeltaF/F';


%% Normalize to baseline mean and std
bl_data = [event_matrix(:, 1:(bl_length)*nsr-1)];
blm = mean(bl_data, 2, 'omitnan'); %Baseline mean pr event
blstd = std(bl_data,[], 2,'omitnan'); %Baseline std pr event
evtmatr_blcor = (event_matrix-blm)./blstd; %Z-score compared to baseline mean and std

tvec = time_event_matrix(1:size(event_matrix,2)); %time vec
zs_std = std(evtmatr_blcor,'omitnan'); %Mean
errvec = std(evtmatr_blcor,[],1,'omitnan')./sqrt(size(evtmatr_blcor, 1)); %Standard error of the mean
mvec = mean(evtmatr_blcor,'omitnan'); %Mean

uerr = mvec+errvec; %upper errorbar
lerr = mvec-errvec; %lower errorbar

figure(4)
h1 = subplot(211);
orgSize1 = get(gca, 'Position');
p2 = plot(tvec, evtmatr_blcor','Color',[0.7 0.7 0.7]);
hold on
p1 = patch([tvec, fliplr(tvec)], [lerr, fliplr(uerr)], lcmap(1, :), 'edgecolor','none');
yline(0,':')
p3 = xline(0, '--');
xline(event_lgth*nT, '--')
p4 = plot(tvec, mvec, 'color', 'k','linewidth', 3);
xlim([tvec(1),  tvec(end)])
title(['Z-score calculated from baseline, first ',num2str(bl_length),' s'])
xlabel('Time (s)')
ylabel('z-score')
legend([p3, p2(1), p4, p1], {'TTL on/offset','Trials', 'Mean response', 'SEM'},'location','northwest')

h2 = subplot(212);
imagesc(tvec, 1, evtmatr_blcor)
xlabel('Time (s)')
set(gca,'YTick', 1:size(evtmatr_blcor,2))
ylabel('Trial #')
colormap(cmap)
orgSize2 = get(gca, 'Position');
cb = colorbar;
set(h1, 'Position', orgSize1)
set(h2,'Position', orgSize2)
cb.Label.String = 'z-score';


%% Data output
sz = [size(event_matrix,2),4+size(event_matrix, 1)];
varTypes = repmat("double", [sz(2),1])';
tnames = [];
for ival = 1:size(event_matrix, 1);
    tmp = strcat("Trial_",num2str(ival));
    tnames = [tnames,tmp];
end
varNames = ["Time (s)", tnames, "Mean response", "std", "SEM"];

dff_table = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
zs_table = dff_table;

dff_table.("Time (s)") = tvec';
dff_table(:,2:sz(2)-3) = array2table(event_matrix');
dff_table.("Mean response") = dff_mvec';
dff_table.std = dff_std';
dff_table.SEM = dff_sem';

zs_table.("Time (s)") = tvec';
zs_table(:,2:sz(2)-3) = array2table(evtmatr_blcor');
zs_table.("Mean response") = mvec';
zs_table.std = zs_std';
zs_table.SEM = errvec';

%%
writetable(dff_table,['Output data/DFF_',filename,'.xlsx']);
writetable(zs_table,['Output data/ZScore_',filename,'.xlsx']);