%% PLOT RESULTS OF "OPTICA" STUDY
% olaf.dimigen@hu-berlin.de, 2018
clear, rng('shuffle')
close all

% restoredefaultpath
addpath('M:/Dropbox/_subfunc_master') % functions: stderr und ci ...
addpath M:/Dropbox/eeglab14_1_2b; eeglab; close
addpath M:/Dropbox/_subfunc_master/export_fig_2018 % export_fig

LC        = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3   3.5 4 5 7.5 10 12.5 15 20 25 30];
EXPORTFIG  = 1 % export figures as .eps?
NCL        = 0 % number of contour lines

%% time and axis limits for sac-locked (SRPs)
TIME1   = -200:2:598;
XMIN1   =  -200;
XMAX1   =  500;
YMIN1c1 =  -11;
YMAX1c1 =   11;

%% plot settings for topos
PRAD       = 0.80 % max(max(chanlocs.radius),0.5) is .75 for these experiments
IRAD       = 0.80
HRAD       = 'rim'
PLOTELECS  = 'off' %  'on','off','labels','numbers','ptslabels','ptsnumbers'
XLIM       = [-0.6 0.6] % so nose is not cut off
YLIM       = XLIM;
SCREENSIZE = get( 0, 'Screensize' ); % values

%% load input data
load Z:/OPTICA/results/metrics.mat metric*
load Z:/OPTICA/results/opticaresults.mat % used only for trial count information
load Z:/OPTICA/results/chanlocs.mat chanlocs_reading chanlocs_scenes
load Z:/OPTICA/results/results_for_plotting.mat SCENE* READ*

%% load time-frequency analysis (TFA) results
load('Z:\OPTICA\scenes\results_TFA\tfa_46chans.mat')

%% elecs: scenes
elec_L  = [01 03 05 08 11 12 16 17 18 22 23 27 28 31 32 36 37 41 46];
elec_C  = [06 09 13 24 33 38 42 44];
elec_R  = [02 04 07 10 14 15 19 20 21 25 26 29 30 34 35 39 40 43 45];

% #########################################################################
% COMPARE 3 CORRECTIONS
% #########################################################################

LWC = 1.2 % line width for central electrodes
CHANS_EEG = 1:46
% % time range and color limit for (CR) topographies
SMP1   = 106 % = from   10 ms after sac onset
SMP2   = 200 % 200 % = until 200 ms after sac onset
CRLIM  =   8 % µV for topoplots
CRLIM2 =   8 % uncorrected
LC1    =   6
LC2    =   8 

fx = figure('name','The bad, the ugly, and the good')

%% TOPOS ------------------------------------------------------------------
% TOPOS: SP (-2 ms)
SMP1 = 100; SMP2 = 100;
subplot(4,9,1); title('-2 ms')
d = mean(SCENE_SACC_ORIG(CHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM2 CRLIM2],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,4); title('-2 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,1,LC1,1),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,7); title('-2 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,2,LC2,2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

% TOPOS: CR (50 ms)    
SMP1 = 130; SMP2 = 130;
subplot(4,9,2); title('60 ms')
d = mean(SCENE_SACC_ORIG(CHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM2 CRLIM2],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,5); title('60 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,1,LC1,1),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,8); title('60 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,2,LC2,2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

% TOPOS: lambda response (120 ms)
SMP1 = 100+61; SMP2 = 100+61;
subplot(4,9,3); title('120 ms')
d = mean(SCENE_SACC_ORIG(CHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM2 CRLIM2],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,6); title('120 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,1,LC1,1),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,9,9); title('120 ms')
d = mean(SCENE_SACC_CORR(CHANS_EEG,SMP1:SMP2,2,LC2,2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
originalSize = get(gca, 'Position'); colorbar; set(gca, 'Position', originalSize);

%% SACCADE-ERPs -----------------------------------------------------------

% BAD
subplot(4,9,10:12); hold on; title('uncorrected')
% plot([0 0],[-80 80],'k:')
plot([  -2 -2],[-15 15],'k')
plot([  60 60],[-15 15],'k')
plot([120 120],[-15 15],'k')
plot(TIME1,SCENE_SACC_ORIG(elec_L,:,1),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_ORIG(elec_R,:,1),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_ORIG(elec_C,:,1),'k','linewidth',LWC)
plot(TIME1,SCENE_SACC_ORIG(42,:,1),'g','linewidth',3)

xlim([XMIN1 XMAX1]), set(gca,'ticklength',1.2*get(gca,'ticklength'))
ylim([-90 90])
set(gca,'Ytick',-80:40:80)
ylabel({'Saccade-related potential [µV]'});
xlabel('Time after saccade (ms)')
LC1 = 6 % 1.0 Hz
LC2 = 8 % 2.5 Hz

% UGLY
subplot(4,9,13:15); hold on; title('Typical')
% plot([0 0],[-80 80],'k:')
plot([  -2 -2],[-15 15],'k')
plot([  60 60],[-15 15],'k')
plot([120 120],[-15 15],'k')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,1,LC1,2),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,1,LC1,1),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,1,LC1,1),'k','linewidth',LWC)
plot(TIME1,SCENE_SACC_CORR(42,:,1,LC1,1),'g','linewidth',3)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',1.2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1])
set(gca,'Ytick',-50:10:50)
xlabel('Time after saccade (ms)')

% GOOD
subplot(4,9,16:18); hold on; title('Optimized')
% plot([0 0],[-80 80],'k:')
plot([  -2 -2],[-15 15],'k')
plot([  60 60],[-15 15],'k')
plot([120 120],[-15 15],'k')
h1 = plot(TIME1,SCENE_SACC_CORR(elec_L,:,2,LC2,2),'color',[0.3 0.3 1.0])
h2 = plot(TIME1,SCENE_SACC_CORR(elec_R,:,2,LC2,2),'color',[1.0 0.3 0.3])
h3 = plot(TIME1,SCENE_SACC_CORR(elec_C,:,2,LC2,2),'k','linewidth',LWC)
h0 = plot(TIME1,SCENE_SACC_CORR(42,    :,2,LC2,2),'g','linewidth',3)

xlim([XMIN1 XMAX1]), set(gca,'ticklength',1.2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1]), % set(gca,'YColor','none')
set(gca,'Ytick',-50:10:50)
xlabel('Time after saccade (ms)')
ll = legend([h0 h1(1) h2(1) h3(1)],{'Oz','left hemisph.','right hemisph.','midline'})
set(ll,'box','off')


%% TFA --------------------------------------------------------------------
freqs = tfa.subj(1).lc(8).hc(2).ow(2).freqs; % get freq levels once
times = tfa.subj(1).lc(8).hc(2).ow(2).times; % get time points once

%1:4 EOG
% 42 Oz
% 9  Afz
% 6  Fpz

CHAN = 6 %[42] % do this plot for Fpz and for electrode Oz

% get subject-level TFA results for raw data, "typical" ICA and "optimal" ICA
for s = 1:12
   ersp1(:,:,:,s) = tfa.subj(s).ersp_raw;  % raw
   ersp2(:,:,:,s) = tfa.subj(s).lc(6).hc(1).ow(1).ersp_corr; % typical
   ersp3(:,:,:,s) = tfa.subj(s).lc(8).hc(2).ow(2).ersp_corr; % near-optimal    
end

% average across subjects (and channels, if multiple)
ga_ersp1 = squeeze(mean(mean(ersp1(CHAN,:,:,:),4),1));
ga_ersp2 = squeeze(mean(mean(ersp2(CHAN,:,:,:),4),1));
ga_ersp3 = squeeze(mean(mean(ersp3(CHAN,:,:,:),4),1));

CLIMS1 = [-10 10];
CLIMS2 = [-2.5 2.5];
% CLIMS2 = [-5 5];
fx_fpz = figure

% figure
subplot(4,9,19:21); % title('uncorrected')
imagesc(times,freqs,ga_ersp1); set(gca,'YDir','normal'); colormap jet; caxis(CLIMS2); xlim([-200 600]);
set(gca,'Ytick',[10:10:50]); set(gca,'XTick',[-100:100:400]);
xlim([XMIN1 XMAX1])
xlabel('Time after saccade (ms)'); ylabel({'Frequency [Hz]'});

subplot(4,9,22:24); % title('Typical ICA')
imagesc(times,freqs,ga_ersp2); set(gca,'YDir','normal'); colormap jet; caxis(CLIMS2); xlim([-200 600]);
set(gca,'Ytick',[10:10:50]); set(gca,'XTick',[-100:100:400]);
xlim([XMIN1 XMAX1])
xlabel('Time after saccade (ms)'); %ylabel({'Frequency [Hz]'});

subplot(4,9,25:27); % title('Optimized ICA')
imagesc(times,freqs,ga_ersp3); set(gca,'YDir','normal'); colormap jet; caxis(CLIMS2); xlim([-200 600]);
set(gca,'Ytick',[10:10:50]); set(gca,'XTick',[-100:100:400]);
xlim([XMIN1 XMAX1])
xlabel('Time after saccade (ms)'); %ylabel({'Frequency [Hz]'});
originalSize = get(gca, 'Position'); hcb = colorbar; set(gca, 'Position', originalSize);
ylabel(hcb,'ERSP (dB)')

% uncorrected minus corrected
subplot(4,9,28:30)
imagesc(times,freqs,ga_ersp1 - ga_ersp3); set(gca,'YDir','normal'); colormap jet; caxis(CLIMS2); xlim([-200 600]);
xlim([XMIN1 XMAX1])
xlabel('Time after saccade (ms)'); ylabel({'Frequency [Hz]'});


%% export figure
if EXPORTFIG
    set(fx_fpz,'Position',SCREENSIZE)
    export_fig(fx_fpz,'figs_2019_R1/figure_COMPARE3.pdf','-painters','-transparent') %  %  '-eps',
end