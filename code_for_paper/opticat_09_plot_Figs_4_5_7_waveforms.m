%% PLOT RESULTS OF "OPTICA" STUDY
% olaf.dimigen@hu-berlin.de, 2018
clear, rng('shuffle')
close all

restoredefaultpath
addpath('M:/Dropbox/_subfunc_master') % functions: stderr und ci ...
addpath M:/Dropbox/eeglab14_1_2b; eeglab; close
addpath M:/Dropbox/_subfunc_master/export_fig_2018 % export_fig

LC        = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3   3.5 4 5 7.5 10 12.5 15 20 25 30];
LCx       = {'No' '0.1 Hz' '0.25 Hz' '0.5 Hz' '0.75 Hz' '1 Hz' '1.5 Hz' '2 Hz' '2.5 Hz' '3 Hz'   '3.5 Hz' '4 Hz' '5 Hz' '7.5 Hz' '10 Hz' '12.5 Hz' '15 Hz' '20 Hz' '25 Hz' '30 Hz'};
PLOTFREQS = [1 3 4 6 8 10 13 17]
HC        = [40 100];
THR       = [0.5:0.1:1.6]; % thresholds, 1.6 = outlier-detection

SUBJECTS   = 1:12
NCHANS_EEG = 46
EXPORTFIG  = 0 % export figures as .eps?

NCL        = 0 % number of contour lines

%% time and axis limits for sac-locked (SRPs)
TIME1   = -200:2:598;
XMIN1   = -100;
XMAX1   =  250;
YMIN1uc =  -80;
YMAX1uc =   80;
YMIN1c1 =  -50;
YMAX1c1 =   50;
YMIN1c2  = -20;
YMAX1c2  =  20;
YMIN1c3  = -20;
YMAX1c3  =  20;
% and for stim-locked (ERPs)
TIME2   = -100:2:298;
XMIN2   = -100
XMAX2   =  200 % saccade-free interval

%% plot settings for topos
PRAD       = 0.80 % max(max(chanlocs.radius),0.5) is .75 for these experiments
IRAD       = 0.80
HRAD       = 'rim'
PLOTELECS  = 'off' %  'on','off','labels','numbers','ptslabels','ptsnumbers'
XLIM       = [-0.6 0.6] % so nose is not cut off
YLIM       = XLIM;
SCREENSIZE = get( 0, 'Screensize' ); % values

%% load input data
load Y:/OPTICA/results/metrics.mat metric*
load Y:/OPTICA/results/opticaresults.mat % used only for trial count information
load Y:/OPTICA/results/chanlocs.mat chanlocs_reading chanlocs_scenes
load Y:/OPTICA/results/results_for_plotting.mat SCENE* READ*


% #########################################################################
% MONSTERFIGURE: SACCADE-LOCKED WAVEFORMS
% #########################################################################

% -- SUMMARY: Aggregated results at threshold of 1.1 (filters x overweighting)
THRESHOLD2PLOT = 7
LW             = 1.0; % lines
LW2            = 0.7; % error bars

% linecolors1
linecolors1 = [0         0    0.5625;
         0         0    0.8125;
         0    0.0625    1.0000;
         0    0.3125    1.0000;
         0    0.5625    1.0000;
         0    0.8125    1.0000;
    0.0625    1.0000    0.9375;
    0.3125    1.0000    0.6875;
    0.5625    1.0000    0.4375;
    0.8125    1.0000    0.1875;
    1.0000    0.9375         0;
    1.0000    0.6875         0;
    1.0000    0.4375         0;
    1.0000    0.1875         0;
    0.9375         0         0;
    0.6875         0         0];

LWC = 0.8 % line width for central electrodes

% time range and color limit for (CR) topographies
SMP1  = 106 % = from   10 ms after sac onset
SMP2  = 200 %200 % = until 200 ms after sac onset
CRLIM =   4 % µV for topoplots

% which version of data to plot
hc    = 2 % 40 or 100 Hz low-pass filter
ow    = 2 % basic/overweighted

fx = figure('name','Saccade-locked ERP')

% elecs: scenes
elec_L  = [01 03 05 08 11 12 16 17 18 22 23 27 28 31 32 36 37 41 46];
elec_C  = [06 09 13 24 33 38 42 44];
elec_R  = [02 04 07 10 14 15 19 20 21 25 26 29 30 34 35 39 40 43 45];

% WAVEFORMS: SCENES
subplot(4,10,1); hold on; title('uncorrected')
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_ORIG(elec_L,:,1),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_ORIG(elec_R,:,1),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_ORIG(elec_C,:,1),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1uc YMAX1uc])
set(gca,'Ytick',-80:40:80)
ylabel({'Uncorrected [µV]'});

subplot(4,10,2); hold on; title(LCx(PLOTFREQS(1)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(1),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(1),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(1),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1])
ylabel({'Corrected [µV]'});

subplot(4,10,3); hold on; title(LCx(PLOTFREQS(2)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(2),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(2),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(2),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1]), set(gca,'YColor','none')
set(gca,'YColor','none'); % remove Y-axis

subplot(4,10,4); hold on; title(LCx(PLOTFREQS(3)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(3),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(3),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(3),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c2 YMAX1c2]), %set(gca,'YColor','none')

subplot(4,10,5); hold on; title(LCx(PLOTFREQS(4)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(4),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(4),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(4),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c2 YMAX1c2]), set(gca,'YColor','none')

subplot(4,10,6); hold on; title(LCx(PLOTFREQS(5)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(5),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(5),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(5),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,7); hold on; title(LCx(PLOTFREQS(6)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(6),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(6),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(6),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,8); hold on; title(LCx(PLOTFREQS(7)))
plot([0 0],[-80 80],'k:')
plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(7),ow),'color',[0.3 0.3 1.0])
plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(7),ow),'color',[1.0 0.3 0.3])
plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(7),ow),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,9); hold on; title(LCx(PLOTFREQS(8)))
plot([0 0],[-80 80],'k:')
h1 = plot(TIME1,SCENE_SACC_CORR(elec_L,:,hc,PLOTFREQS(8),ow),'color',[0.3 0.3 1.0]);
h2 = plot(TIME1,SCENE_SACC_CORR(elec_R,:,hc,PLOTFREQS(8),ow),'color',[1.0 0.3 0.3]);
h3 = plot(TIME1,SCENE_SACC_CORR(elec_C,:,hc,PLOTFREQS(8),ow),'k','linewidth',LWC);
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')
l = legend([h1(1) h2(1) h3(1)],{'L hemisph.','R hemisph.','Midline'},'box','off')

subplot(4,10,10); hold on; title('MSEC')
plot(TIME1,mean(SCENE_SACC_CORR_MSEC(elec_L,:,1:12),3),'color',[0.3 0.3 1.0])
plot(TIME1,mean(SCENE_SACC_CORR_MSEC(elec_R,:,1:12),3),'color',[1.0 0.3 0.3])
plot(TIME1,mean(SCENE_SACC_CORR_MSEC(elec_C,:,1:12),3),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')


% TOPOS: SCENES
subplot(4,10,11);
d = mean(SCENE_SACC_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
%originalSize = get(gca, 'Position'); colorbar; set(gca, 'Position', originalSize);

subplot(4,10,12);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(1),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,13);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(2),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,14);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(3),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,15);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(4),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,16);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(5),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,17);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(6),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,18);
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(7),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,19); 
d = mean(SCENE_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(8),ow),2);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,20); % MSEC
d = mean(mean(SCENE_SACC_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2,:),2),3);
topoplot(d,chanlocs_scenes,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
originalSize = get(gca, 'Position'); colorbar; set(gca, 'Position', originalSize);


% WAVEFORMS: READING
% elecs: reading
elec_L = [01 03 05 08 11 12 16 17 18 22 23 27 28 32 33 37 38 42 46]
elec_R = [02 04 07 10 14 15 19 20 21 25 26 29 30 31 35 36 40 41 44]
elec_C = [06 09 13 24 34 39 43 45]

subplot(4,10,21); hold on; title('uncorrected')
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_ORIG(elec_L,:,1),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_ORIG(elec_R,:,1),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_ORIG(elec_C,:,1),'k','linewidth',LWC)
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1uc YMAX1uc])
set(gca,'Ytick',-80:40:80)
ylabel({'Uncorrected [µV]'});

subplot(4,10,22); hold on; title(LCx(PLOTFREQS(1)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(1),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(1),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(1),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,1),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1])
ylabel({'Corrected [µV]'});

subplot(4,10,23); hold on; title(LCx(PLOTFREQS(2)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(2),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(2),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(2),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,2),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c1 YMAX1c1]), set(gca,'YColor','none')
set(gca,'YColor','none'); % remove Y-axis

subplot(4,10,24); hold on; title(LCx(PLOTFREQS(3)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(3),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(3),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(3),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,3),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c2 YMAX1c2]), %set(gca,'YColor','none')

subplot(4,10,25); hold on; title(LCx(PLOTFREQS(4)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(4),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(4),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(4),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,4),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c2 YMAX1c2]), set(gca,'YColor','none')

subplot(4,10,26); hold on; title(LCx(PLOTFREQS(5)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(5),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(5),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(5),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,5),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,27); hold on; title(LCx(PLOTFREQS(6)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(6),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(6),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(6),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,6),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,28); hold on; title(LCx(PLOTFREQS(7)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(7),ow),'color',[0.3 0.3 1.0])
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(7),ow),'color',[1.0 0.3 0.3])
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(7),ow),'k','linewidth',LWC)
% plot(TIME1,READ_SACC_CORR(elec_SP,:,7),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,29); hold on; title(LCx(PLOTFREQS(8)))
plot([0 0],[-80 80],'k:')
plot(TIME1,READ_SACC_CORR(elec_L,:,hc,PLOTFREQS(8),ow),'color',[0.3 0.3 1.0]);
plot(TIME1,READ_SACC_CORR(elec_R,:,hc,PLOTFREQS(8),ow),'color',[1.0 0.3 0.3]);
plot(TIME1,READ_SACC_CORR(elec_C,:,hc,PLOTFREQS(8),ow),'k','linewidth',LWC);
%legend([h11(1) h12(1) h13(1)],{'L hemisph.','R hemisph.','midline'},'box','off')
% plot(TIME1,READ_SACC_CORR(elec_SP,:,8),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

subplot(4,10,30); hold on; title('MSEC')
plot(TIME1,mean(READ_SACC_CORR_MSEC(elec_L,:,1:12),3),'color',[0.3 0.3 1.0])
plot(TIME1,mean(READ_SACC_CORR_MSEC(elec_R,:,1:12),3),'color',[1.0 0.3 0.3])
plot(TIME1,mean(READ_SACC_CORR_MSEC(elec_C,:,1:12),3),'k','linewidth',LWC)
% plot(TIME1,mean(READ_SACC_CORR_MSEC(elec_SP,:,1:12),3),'color',[1.0 0.3 0.3])
xlim([XMIN1 XMAX1]), set(gca,'ticklength',2*get(gca,'ticklength'))
ylim([YMIN1c3 YMAX1c3]), set(gca,'YColor','none')

% TOPOS: READING
subplot(4,10,31); 
d = mean(READ_SACC_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
%originalSize = get(gca, 'Position'); colorbar; set(gca, 'Position', originalSize);

subplot(4,10,32);
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(1),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,33); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(2),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,34); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(3),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,35);
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(4),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,36); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(5),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,37); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(6),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,38); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(7),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,39); 
d = mean(READ_SACC_CORR(1:NCHANS_EEG,SMP1:SMP2,hc,PLOTFREQS(8),ow),2);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,40); 
d = mean(mean(READ_SACC_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2,:),2),3);
topoplot(d,chanlocs_reading,'maplimits',[-CRLIM CRLIM],'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
originalSize = get(gca, 'Position'); colorbar; set(gca, 'Position', originalSize);
axis([XLIM YLIM])
    
% export figure
if EXPORTFIG
    set(gcf,'Position',SCREENSIZE)
    export_fig(gcf,'figs_2018/figure_waves_saccadelocked_new2.eps','-eps','-transparent','-painters')
end

% #########################################################################
% SPIKE POTENTIAL WAVEFORMS
% #########################################################################
% PLOTFREQS = [1 4 6 8 10 14 17 20]
PLOTFREQS_SP = [1 4 6 8 10]

LW = 1.1;

% from colormap "hot"
linecolors = [0.0416  0  0;
                0.375  0  0;
                0.708333  0  0;
                1        0.0416  0;
                1        0.375  0;
                1        0.7083  0;
                1  1     0.0625;
                1  1     0.5625];


%% ALL-IN-ONE PLOT SPIKE POTENTIAL MAARTEN-STYLE
fig_sp = figure; 

CLUSTER_A = [1:4]; % electrode cluster 1
CLUSTER_B = 38;    % electrode cluster 2 
YLIMITZ1 = [-32 15];
YLIMITZ2 = [-24 15];
XLIMITZ  = [-35 30];

% SCENES
subplot(2,5,1); hold on; title('Uncorrected')
plot([0 0],[-80 80],'k:')
plot(TIME1,mean(SCENE_SACC_ORIG(CLUSTER_A,:),1)-mean(SCENE_SACC_ORIG(CLUSTER_B,:),1),'color',[0.0 0.0 0.0],'linewidth',LW)
plot(TIME1,mean(mean(SCENE_SACC_CORR_MSEC(CLUSTER_A,:,1:12),1),3)-mean(mean(SCENE_SACC_CORR_MSEC(CLUSTER_B,:,1:12),1),3),'color',[0.0 0.0 1.0],'linewidth',LW,'linestyle','--')

ylim(YLIMITZ1)
xlim(XLIMITZ)
ylabel('Spike potential: Scenes [µV]')
xlabel('Time after saccade [ms]')

subplot(2,5,2); hold on; title('40 Hz | Basic')
plot([0 0],[-80 80],'k:'); hold on
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(SCENE_SACC_CORR(CLUSTER_A,:,1,PLOTFREQS_SP(f),1),1)-mean(SCENE_SACC_CORR(CLUSTER_B,:,1,PLOTFREQS_SP(f),1),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ1)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

subplot(2,5,3); hold on; title('40 Hz | Overweighted')
plot([0 0],[-80 80],'k:'); hold on
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(SCENE_SACC_CORR(CLUSTER_A,:,1,PLOTFREQS_SP(f),2),1)-mean(SCENE_SACC_CORR(CLUSTER_B,:,1,PLOTFREQS_SP(f),2),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ1)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

subplot(2,5,4); hold on; title('100 Hz | Basic')
plot([0 0],[-80 80],'k:'); hold on
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(SCENE_SACC_CORR(CLUSTER_A,:,2,PLOTFREQS_SP(f),1),1)-mean(SCENE_SACC_CORR(CLUSTER_B,:,2,PLOTFREQS_SP(f),1),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ1)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

subplot(2,5,5); hold on; title('100 Hz | Overweighted')
plot([0 0],[-80 80],'k:')
for f = 1:length(PLOTFREQS_SP)
    h(f) = plot(TIME1,mean(SCENE_SACC_CORR(CLUSTER_A,:,2,PLOTFREQS_SP(f),2),1)-mean(SCENE_SACC_CORR(CLUSTER_B,:,2,PLOTFREQS_SP(f),2),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ1)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')
%l = legend(h,{LCx{PLOTFREQS_SP(1)},LCx{PLOTFREQS_SP(2)},LCx{PLOTFREQS_SP(3)},LCx{PLOTFREQS_SP(4)},LCx{PLOTFREQS_SP(5)},LCx{PLOTFREQS_SP(6)},LCx{PLOTFREQS_SP(7)},LCx{PLOTFREQS_SP(8)} })
l = legend(h,{'No',LCx{PLOTFREQS_SP(2)},LCx{PLOTFREQS_SP(3)},LCx{PLOTFREQS_SP(4)},LCx{PLOTFREQS_SP(5)} })
set(l,'Location','South','box','off')

% READING
subplot(2,5,6); hold on; title('Uncorrected')
plot([0 0],[-80 80],'k:')
plot(TIME1,mean(READ_SACC_ORIG(CLUSTER_A,:),1)-mean(READ_SACC_ORIG(CLUSTER_B,:),1),'color',[0.0 0.0 0.0],'linewidth',1.2,'linewidth',LW)
plot(TIME1,mean(mean(READ_SACC_CORR_MSEC(CLUSTER_A,:,1:12),1),3)-mean(mean(READ_SACC_CORR_MSEC(CLUSTER_B,:,1:12),1),3),'color',[0.0 0.0 1.0],'linewidth',LW,'linestyle','--')
ylim(YLIMITZ2)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')
ylabel('Spike potential: Reading [µV]')

subplot(2,5,7); hold on; title('40 Hz | Basic')
plot([0 0],[-80 80],'k:')
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(READ_SACC_CORR(CLUSTER_A,:,1,PLOTFREQS_SP(f),1),1)-mean(READ_SACC_CORR(CLUSTER_B,:,1,PLOTFREQS_SP(f),1),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ2)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')
ylabel('Spike potential [µV]')

subplot(2,5,8); hold on; title('40 Hz | Overweighted')
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(READ_SACC_CORR(CLUSTER_A,:,1,PLOTFREQS_SP(f),2),1)-mean(READ_SACC_CORR(CLUSTER_B,:,1,PLOTFREQS_SP(f),2),1),'color',linecolors(f,:),'linewidth',LW)
end
plot([0 0],[-80 80],'k:')
ylim(YLIMITZ2)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

subplot(2,5,9); hold on; title('100 Hz | Basic')
plot([0 0],[-80 80],'k:')
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(READ_SACC_CORR(CLUSTER_A,:,2,PLOTFREQS_SP(f),1),1)-mean(READ_SACC_CORR(CLUSTER_B,:,2,PLOTFREQS_SP(f),1),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ2)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

subplot(2,5,10); hold on; title('100 Hz | Overweighted')
plot([0 0],[-80 80],'k:')
for f = 1:length(PLOTFREQS_SP)
    plot(TIME1,mean(READ_SACC_CORR(CLUSTER_A,:,2,PLOTFREQS_SP(f),2),1)-mean(READ_SACC_CORR(CLUSTER_B,:,2,PLOTFREQS_SP(f),2),1),'color',linecolors(f,:),'linewidth',LW)
end
ylim(YLIMITZ2)
xlim(XLIMITZ)
xlabel('Time after saccade [ms]')

if EXPORTFIG
    % export figure
    set(fig_sp,'Position',[50 50 1100 650])
    export_fig(fig_sp,'figs_2018/figure_waves_spikepotential','-eps','-transparent','-painters')
end


% #########################################################################
% STIM-ERP and EYE-TRACK
% #########################################################################
ow = 1

fig_erp = figure

XMAX2        = 200
MYCHAN_READ  = 41 % "PO10" for reading
MYCHAN_SCENE = 42 % "Oz" for scenes
% MYCHAN       =  1

% scenes: eye-track & ERP differences
ET_SCENE = rmbase(SCENE_STIM_ORIG(NCHANS_EEG+1:NCHANS_EEG+4,:,1,ow),89,[1:50]).* 0.036;  % 89 ? ###
ET_READ  = rmbase(READ_STIM_ORIG(NCHANS_EEG+1:NCHANS_EEG+4,:,1,ow),89,[1:50]).* 0.036;  % 89 ? ###
ET_SCENE = SCENE_STIM_ORIG(NCHANS_EEG+1:NCHANS_EEG+4,:,1,ow).* 0.036;
ET_READ  = READ_STIM_ORIG(NCHANS_EEG+1:NCHANS_EEG+4,:,1,ow) .* 0.036;

subplot(2,2,[1]); hold on; title('Scenes: Image onset')
plot([0 0],[-15 22],'k')
plot(TIME2,SCENE_STIM_ORIG(1:NCHANS_EEG,:,1,ow),'k')
plot(TIME2,SCENE_STIM_ORIG(MYCHAN_SCENE,:,1,ow),'r','linewidth',1.8)
%plot(TIME2,SCENE_STIM_ORIG([2 7 12 17 22 27 32 37 42],:,1,ow),'g','linewidth',1.8)
% plot(TIME2,SCENE_STIM_ORIG([42],:,1,ow),'b','linewidth',1.8)
% plot(TIME2,SCENE_STIM_ORIG([2 7 12 17 22],:,1,ow),'g','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 22])
ylabel({'Stimulus-ERP [µV]'})
xlabel('Time after stimulus onset (ms)')

subplot(2,2,[3]); hold on; title('Eye-Track')
plot([0 0],[-2 2],'k')
h1 = plot(TIME2,ET_SCENE([1],:),'b');
h2 = plot(TIME2,ET_SCENE([2],:),'b--');
h3 = plot(TIME2,ET_SCENE([3],:),'g');
h4 = plot(TIME2,ET_SCENE([4],:),'g--');
% h1 = plot(TIME2,ET_SCENE([1 3],:),'b');
% h2 = plot(TIME2,ET_SCENE([2 4],:),'c');
legend([h1 h2 h3 h4],{'left eye horiz.','left eye vertic.','right eye horiz.','right eye vertic.'})
xlim([XMIN2 XMAX2])
ylim([-0.3 0.3])
ylabel('Gaze position (deg)')
xlabel('Time after stimulus onset (ms)')

subplot(2,2,[2]); hold on; title('Reading: Sentence onset')
plot([0 0],[-15 22],'k')
plot(TIME2,READ_STIM_ORIG(1:NCHANS_EEG,:,1,ow),'k')
plot(TIME2,READ_STIM_ORIG(MYCHAN_READ,:,1,ow),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 22])
ylabel({'Stimulus-ERP [µV]'})
xlabel('Time after stimulus onset (ms)')

subplot(2,2,[4]); hold on; title('Eye-Track')
plot([0 0],[-2 2],'k')
h1 = plot(TIME2,ET_READ([1],:),'b');
h2 = plot(TIME2,ET_READ([2],:),'b--');
h3 = plot(TIME2,ET_READ([3],:),'g');
h4 = plot(TIME2,ET_READ([4],:),'g--');
% h1 = plot(TIME2,ET_READ([1 3],:),'b');
% h2 = plot(TIME2,ET_READ([2 4],:),'c');
% legend([h1(1) h2(1)],{'hor. gaze','vert. gaze'})
xlim([XMIN2 XMAX2])
ylim([-0.3 0.3])
ylabel('gaze position (deg)')
xlabel('Time after stimulus onset (ms)')
% fig_erp


%% get single-trial gaze data
load Y:\OPTICA\results\singletrialgaze.mat 

LX1    = rmbase(singletrialgazeLX1.*0.036,150,1:150); % single-trial
SD_LX1 = std(LX1,0,1); % SD of single-trials
LX1m    = mean(LX1,1); % mean

LY1    = rmbase(singletrialgazeLY1.*0.036,150,1:150);
SD_LY1 = std(LY1,0,1);
LY1m    = mean(LY1,1);

figure; hold on; title('single-trial gaze data (baselined) plus SD')
plot(LX1m,'r');
plot(LX1m+SD_LX1,'c:')
plot(LX1m-SD_LX1,'c:')


% #########################################################################
%                         OVERCORRECTION WAVEFORMS                            
% #########################################################################

HC        = 2
OW        = 2
PLOTFREQS = [1 3 4 6 8 10 13 17 18]
SMP1      = 50     % same window as in select_IC (technically, should be 51)
SMP2      = 50+100 % same window as in select_IC 
MLIM1     = -1.5
MLIM2     = +1.5

f5 = figure('name','scenes: stim-locked')

subplot(4,10,1); hold on; title('Stim-ERP')
plot(TIME2,SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-20 20],'k')
ylabel('Stimulus-ERP [µV] ')

subplot(4,10,2); hold on; title(LCx(PLOTFREQS(1)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,1)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(1),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel('ERP distortion [µV] ')

subplot(4,10,3); hold on; title(LCx(PLOTFREQS(2)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(2),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(2),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,4); hold on; title(LCx(PLOTFREQS(3)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(3),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(3),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,5); hold on; title(LCx(PLOTFREQS(4)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(4),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(4),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,6); hold on; title(LCx(PLOTFREQS(5)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(5),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(5),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,7); hold on; title(LCx(PLOTFREQS(6)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(6),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(6),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,8); hold on; title(LCx(PLOTFREQS(7)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(7),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(7),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,9); hold on; title(LCx(PLOTFREQS(8)))
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(8),OW)-SCENE_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN_SCENE,:,HC,PLOTFREQS(8),OW)-SCENE_STIM_ORIG(MYCHAN_SCENE,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,10); hold on; title('MSEC')
plot(TIME2,SCENE_STIM_ORIG(1:NCHANS_EEG,:)-mean(SCENE_STIM_CORR_MSEC(1:NCHANS_EEG,:,1:12),3),'k')
plot(TIME2,SCENE_STIM_ORIG(MYCHAN_SCENE,:)-mean(SCENE_STIM_CORR_MSEC(MYCHAN_SCENE,:,1:12),3),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

%% SCENES: STIM: TOPOS

% absolute
DX = mean(SCENE_STIM_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);

subplot(4,10,11);  hold on
% topoplot(DX,chanlocs_scenes,'maplimits',[-15 15], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
% axis([XLIM YLIM])
plot([0 0],[-2 2],'k')
% h1 = plot(TIME2,ET_SCENE([1],:),'k');
% h2 = plot(TIME2,ET_SCENE([2],:),'k--');
% h3 = plot(TIME2,ET_SCENE([3],:),'b');
% h4 = plot(TIME2,ET_SCENE([4],:),'b--');
h1 = plot(TIME2,mean(ET_SCENE([1 3],:),1),'k');
h2 = plot(TIME2,mean(ET_SCENE([2 4],:),1),'b');
xlim([XMIN2 XMAX2])
ylim([-0.3 0.3])
ylabel('Gaze position (deg)')
xlabel('Time after stimulus onset (ms)')
%legend([h1 h2 h3 h4],{'Left eye horiz.','Left eye vertic.','Right eye horiz.','Right eye vertic.'})


% differences!
subplot(4,10,12);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(1),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,13);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(2),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,14);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(3),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,15);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(4),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,16);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(5),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,17);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(6),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,18);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(7),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,19);
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(8),OW),2);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
orgSize = get(gca, 'Position'); set(gca, 'Position', orgSize);
axis([XLIM YLIM])

subplot(4,10,20); title('MSEC');  hold on;
d = mean(mean(SCENE_STIM_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2,:),2),3);
topoplot(d-DX,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);
axis([XLIM YLIM])



%% reading
subplot(4,10,21); hold on; title('Stim-ERP')
plot(TIME2,READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel('Stimulu-ERP [µV] ')

subplot(4,10,22); hold on; title(LCx(PLOTFREQS(1)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(1),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(1),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel('ERP distortion [µV] ')

subplot(4,10,23); hold on; title(LCx(PLOTFREQS(2)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(2),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(2),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,24); hold on; title(LCx(PLOTFREQS(3)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(3),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(3),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,25); hold on; title(LCx(PLOTFREQS(4)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(4),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(4),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,26); hold on; title(LCx(PLOTFREQS(5)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(5),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(5),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,27); hold on; title(LCx(PLOTFREQS(6)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(6),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(6),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,28); hold on; title(LCx(PLOTFREQS(7)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(7),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(7),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,29); hold on; title(LCx(PLOTFREQS(8)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(8),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(8),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,30); hold on; title(LCx(PLOTFREQS(9)))
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(9),OW)-READ_STIM_ORIG(1:NCHANS_EEG,:),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN_READ,:,HC,PLOTFREQS(9),OW)-READ_STIM_ORIG(MYCHAN_READ,:),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,30); hold on; title('MSEC')
plot(TIME2,READ_STIM_ORIG(1:NCHANS_EEG,:)-mean(READ_STIM_CORR_MSEC(1:NCHANS_EEG,:,1:12),3),'k')
plot(TIME2,READ_STIM_ORIG(MYCHAN_READ,:)-mean(READ_STIM_CORR_MSEC(MYCHAN_READ,:,1:12),3),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

% absolute
DX = mean(READ_STIM_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);

subplot(4,10,31); hold on
% topoplot(DX,chanlocs_reading,'maplimits',[-5 5], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
% axis([XLIM YLIM])
plot([0 0],[-2 2],'k')
% h1 = plot(TIME2,ET_READ([1],:),'k');
% h2 = plot(TIME2,ET_READ([2],:),'k--');
% h3 = plot(TIME2,ET_READ([3],:),'b');
% h4 = plot(TIME2,ET_READ([4],:),'b--');
h1 = plot(TIME2,mean(ET_READ([1 3],:),1),'k');
h2 = plot(TIME2,mean(ET_READ([2 4],:),1),'b-');
legend([h1(1) h2(1)],{'Mean hor. gaze','Mean vert. gaze'})
xlim([XMIN2 XMAX2])
ylim([-0.3 0.3])
ylabel('gaze position (deg)')
xlabel('Time after stimulus onset (ms)')
% fig_erp

% differences!
subplot(4,10,32)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(1),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,33)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(2),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,34)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(3),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,35)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(4),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,36)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(5),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,37)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(6),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,38)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(7),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,39)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(8),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); set(gca, 'Position', orgSize);

subplot(4,10,40)
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(9),OW),2);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

subplot(4,10,40)
d = mean(mean(READ_STIM_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2,:),2),3);
topoplot(d-DX,chanlocs_reading,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

% % title for figure
% [axt,ht]=suplabel('Reading: Stimulus-ERP distortion','t');
% set(ht,'FontSize',18)
% [axty,hty]=suplabel('Stim-ERP  /  Stim-ERP distortion','y');
% set(hty,'FontSize',18)
% % title for figure
[axty,hty]=suplabel('Overcorrection','t');
set(hty,'FontSize',18)

if EXPORTFIG
    % export figure
    set(f5,'Position',SCREENSIZE)
    export_fig(f5,'figs_2018/bigfig_overcorrection.eps','-eps','-transparent','-painters')
end





% #########################################################################1
% ABSOLUTE STIM-ERP / OVERCORRECTION
f8 = figure('name','xxx')


subplot(4,10,1); hold on; title('Eye-Track ')
h1 = plot(TIME2,ET([1 3],:),'b');
h2 = plot(TIME2,ET([2 4],:),'c');
legend([h1(1) h2(1)],{'hor. gaze','vert. gaze'})
xlim([XMIN2 XMAX2])
ylim([-1.0 1.0])
plot([0 0],[-4 4],'k')
ylabel('Eye-Position (degree)')

subplot(4,10,2); hold on; title('uncorrected')
plot(TIME2,SCENE_STIM_ORIG(1:NCHANS_EEG,:,HC,PLOTFREQS(1),OW),'k')
plot(TIME2,SCENE_STIM_ORIG(MYCHAN,:,HC,PLOTFREQS(1),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel({'Stimulus-ERP','[µV]'})
xlabel('Time after stimulus (ms)')

subplot(4,10,3); hold on; title('No')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,1),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(1),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel('ERP: ICA-corr. [µV] ')

subplot(4,10,4); hold on; title('0.1 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,2),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(2),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,5); hold on; title('0.5 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,3),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(3),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,6); hold on; title('1 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,4),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(4),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,7); hold on; title('2 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,5),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(5),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,8); hold on; title('3 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,6),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(6),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,9); hold on; title('5 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,7),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(7),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,10); hold on; title('10 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,8),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(8),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,11); hold on; title('20 Hz')
plot(TIME2,SCENE_STIM_CORR(1:NCHANS_EEG,:,9),OW),'k')
plot(TIME2,SCENE_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(9),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,11); hold on; title('MSEC')
% % % plot(TIME2,SCENE_STIM_CORR_MSEC(1:NCHANS_EEG,:),OW),'k')
% % % plot(TIME2,SCENE_STIM_CORR_MSEC(MYCHAN,:,HC,PLOTFREQS(OW),'r','linewidth',1.8)
% % % xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')


subplot(4,10,21); hold on; title('Eye-Track ')
h1 = plot(TIME2,ET([1 3],:),'b');
h2 = plot(TIME2,ET([2 4],:),'c');
legend([h1(1) h2(1)],{'hor. gaze','vert. gaze'})
xlim([XMIN2 XMAX2])
ylim([-1.0 1.0])
plot([0 0],[-4 4],'k')
ylabel('Eye-Position (degree)')

subplot(4,10,21); hold on; title('uncorrected')
plot(TIME2,READ_STIM_ORIG(1:NCHANS_EEG,:,HC,PLOTFREQS(1),OW),'k')
plot(TIME2,READ_STIM_ORIG(MYCHAN,:,HC,PLOTFREQS(1),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel({'Stimulus-ERP','[µV]'})
xlabel('Time after stimulus (ms)')

subplot(4,10,22); hold on; title('No')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(1),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(1),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')
ylabel('ERP: ICA-corr. [µV] ')

subplot(4,10,23); hold on; title('0.1 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(2),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(2),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,24); hold on; title('0.5 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(3),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(3),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,25); hold on; title('1 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(4),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(4),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,26); hold on; title('2 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(5),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(5),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,27); hold on; title('3 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(6),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(6),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,28); hold on; title('5 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(7),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(7),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,29); hold on; title('10 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(8),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(8),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,30); hold on; title('20 Hz')
plot(TIME2,READ_STIM_CORR(1:NCHANS_EEG,:,HC,PLOTFREQS(9),OW),'k')
plot(TIME2,READ_STIM_CORR(MYCHAN,:,HC,PLOTFREQS(9),OW),'r','linewidth',1.8)
xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')

subplot(4,10,11); hold on; title('MSEC')
% % % plot(TIME2,READ_STIM_CORR_MSEC(1:NCHANS_EEG,:,HC,PLOTFREQS(OW),'k')
% % % plot(TIME2,READ_STIM_CORR_MSEC(MYCHAN,:,HC,PLOTFREQS(OW),'r','linewidth',1.8)
% % % xlim([XMIN2 XMAX2]), ylim([-15 15])
plot([0 0],[-15 15],'k')



topos...

subplot(4,10,11); title('Uncorrected');  hold on;
d = mean(SCENE_STIM_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca,'Position'); colorbar; set(gca,'Position', orgSize);

subplot(4,10,12); title('No');  hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(1),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,13); title('0.1 Hz');  hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(2),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,14); title('0.5 Hz'); hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(3),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,15); title('1 Hz'); hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(4),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,16); title('2 Hz'); hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(5),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,17); title('3 Hz'); hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(6),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,18); title('5 Hz'); hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(7),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,19); title('10 Hz');  hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(8),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);
axis([XLIM YLIM])

subplot(4,10,20); title('20 Hz');  hold on;
d = mean(SCENE_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(9),OW),2);
topoplot(d,chanlocs_scenes,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

% % subplot(4,10,10); title('MSEC');  hold on;
% % d = mean(SCENE_STIM_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(8),2);
% % topoplot(d,chanlocs_scenes,'maplimits',[MLIM1 MLIM2], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
% % orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);





subplot(4,10,1); title('Uncorrected');  hold on;
d = mean(READ_STIM_ORIG(1:NCHANS_EEG,SMP1:SMP2),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

subplot(4,10,2); title('No');  hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(1),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,3); title('0.1 Hz');  hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(2),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,4); title('0.5 Hz'); hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(3),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,5); title('1 Hz'); hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(4),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,6); title('2 Hz'); hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(5),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,7); title('3 Hz'); hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(6),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,8); title('5 Hz'); hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(7),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])

subplot(4,10,9); title('10 Hz');  hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(8),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

subplot(4,10,10); title('20 Hz');  hold on;
d = mean(READ_STIM_CORR(1:NCHANS_EEG,SMP1:SMP2,HC,PLOTFREQS(9),2);
topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
axis([XLIM YLIM])
orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);

subplot(4,10,10); title('MSEC');  hold on;
% d = mean(READ_STIM_CORR_MSEC(1:NCHANS_EEG,SMP1:SMP2),2);
% topoplot(d,chanlocs_reading,'maplimits',[-8 8], 'plotrad',PRAD,'intrad',IRAD,'headrad',HRAD,'electrodes',PLOTELECS,'numcontour',NCL);
% orgSize = get(gca, 'Position'); colorbar; set(gca, 'Position', orgSize);