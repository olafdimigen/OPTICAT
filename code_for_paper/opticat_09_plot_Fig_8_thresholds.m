%% PLOT RESULTS OF "OPTICA" STUDY
% olaf.dimigen@hu-berlin.de, 2018
% updated and fixed issue with a non-intuitive Y-scale of Figure 8A (Nov. 2019)

clear, rng('shuffle')
close all

addpath('M:/Dropbox/_subfunc_master') % functions: stderr und ci ...
addpath M:/Dropbox/_subfunc_master/export_fig_2018 % export_fig
addpath M:/Dropbox/_subfunc_master/ufo/newplot % putnewaxes by Ulrich R.

LC  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30];
HC  = [40 100];
THR = [0.5:0.1:1.6]; % thresholds, 1.6 = outlier-detection

SUBJECTS   = 1:12;
EXPORTFIG  = 1; % export figures as .eps?
SCREENSIZE = [10 10 1000 1000]; % size for figure export (pixel)

%% load input data
load Z:/OPTICA/results/metrics_20180604.mat metric* nbadcomps* exp1_* exp2_* THR* LOW*
load Z:/OPTICA/results/opticaresults_20180604.mat % used only for trial count information
load Z:/OPTICA/results/results_for_plotting_20180604.mat SCENE* READ*
load Z:/OPTICA/results/chanlocs.mat chanlocs_reading chanlocs_scenes

%% Constants
LC_TO_PLOT       = 9          % = 2.0 Hz
LW               = 1.1;
LW2              = 0.8;       % linewidth error bars

XLINE            = 7          % where to put x-line, 7 = 1.1
THRESH2PLOT      = 2:11;      % 0.5:1.5 and "auto"
NCHANS_EEG       = 46;

% X-Ticks & X-Tick labels
t2plot_numbs     = [2:11];
t2plot_names     = {'0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5'};

% plot effects of threshold at what HP-filter/LP-filter/ and overweighting level?
LC_LEVEL         = 8; % 2 Hz
HC_LEVEL         = 2; % 100 Hz
OW_LEVEL         = 2; % yes

f_thresh = figure('name','Threshold choice')

%% compute distribution of variance ratios (for scenes & reading)
all_vr1 = []; all_vr2 = []; all_vr12 = [];

for exp = 1:2
    for s = 1:12
        vr = opticaresults(exp).subj(s).hc(HC_LEVEL).lc(LC_LEVEL).ow(OW_LEVEL).varratio50; % use optimal window: 5 samples (10 ms) before sac on until sac off
        % store distributions
        if exp == 1
            all_vr1 = [all_vr1; vr];
        else
            all_vr2 = [all_vr2; vr];
        end
    end
end

% summarize all values equal to or above this
MAX = 1.5
all_vr1(all_vr1 > MAX) = MAX;
all_vr2(all_vr2 > MAX) = MAX;

MIN = 0.6
all_vr1(all_vr1 < MIN) = MIN;
all_vr2(all_vr2 < MIN) = MIN;

% get histogram data
edges = [0.6:0.1:1.5];
n1 = histc(all_vr1, edges);
n2 = histc(all_vr2, edges);

% divide by 12 (normalize by number of participants, makes more sense --> corrected during Proof stage)
n1 = n1/12;
n2 = n2/12;

%% plot distribution of variance ratios
% scenes
subplot(6,2,1); title('Scenes'); hold on
plot([THRESHOLDS(XLINE) THRESHOLDS(XLINE)],[0 8],'k:');
h1 = bar(edges,n1,'facecolor',[0.1 0.1 0.9]);
ylabel('Number of ICs')
xlabel('Variance ratio threshold used')
set(gca,'XTick',edges);
set(gca,'XTickLabel',{'<= 0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','>= 1.5'});
h1.EdgeColor = 'none'; h1.BarWidth = 0.75;
% set(gca,'XTick',t2plot_numbs);
% set(gca,'XTickLabel',t2plot_names)
ylim([0 25])
xlim([0.5 1.6])
set(gca,'YTick',[0:5:25]);

% reading
subplot(6,2,2); title('Reading'); hold on
plot([THRESHOLDS(XLINE) THRESHOLDS(XLINE)],[0 8],'k:');
h2 = bar(edges,n2,'facecolor',[0.9 0.1 0.1]);
ylabel('Number of ICs')
xlabel('Variance ratio threshold used')
set(gca,'XTick',edges);
set(gca,'XTickLabel',{'<= 0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','>= 1.5'});
h2.EdgeColor = 'none'; h2.BarWidth = 0.75;
% set(gca,'XTick',t2plot_numbs);
% set(gca,'XTickLabel',t2plot_names)
ylim([0 25])
xlim([0.5 1.6])
set(gca,'YTick',[0:5:25]);

%% compute number of rejected ICs per threshold

% scenes
nbadcomps_exp1_lc1_2Hz_hc2_ow1  = squeeze(nbadcomps(1,SUBJECTS,1,LC_LEVEL,1,THRESH2PLOT));
nbadcomps_exp1_lc1_2Hz_hc2_ow2  = squeeze(nbadcomps(1,SUBJECTS,1,LC_LEVEL,2,THRESH2PLOT));
nbadcomps_exp1_lc2_2Hz_hc2_ow1  = squeeze(nbadcomps(1,SUBJECTS,2,LC_LEVEL,1,THRESH2PLOT));
nbadcomps_exp1_lc2_2Hz_hc2_ow2  = squeeze(nbadcomps(1,SUBJECTS,2,LC_LEVEL,2,THRESH2PLOT));

% reading
nbadcomps_exp2_lc1_2Hz_hc2_ow1  = squeeze(nbadcomps(2,SUBJECTS,1,LC_LEVEL,1,THRESH2PLOT));
nbadcomps_exp2_lc1_2Hz_hc2_ow2  = squeeze(nbadcomps(2,SUBJECTS,1,LC_LEVEL,2,THRESH2PLOT));
nbadcomps_exp2_lc2_2Hz_hc2_ow1  = squeeze(nbadcomps(2,SUBJECTS,2,LC_LEVEL,1,THRESH2PLOT));
nbadcomps_exp2_lc2_2Hz_hc2_ow2  = squeeze(nbadcomps(2,SUBJECTS,2,LC_LEVEL,2,THRESH2PLOT));

m1_lc1_ow1  = mean(nbadcomps_exp1_lc1_2Hz_hc2_ow1,1);
m2_lc1_ow1  = mean(nbadcomps_exp2_lc1_2Hz_hc2_ow1,1);
se1_lc1_ow1 =  std(nbadcomps_exp1_lc1_2Hz_hc2_ow1); % alternative stderr()
se2_lc1_ow1 =  std(nbadcomps_exp2_lc1_2Hz_hc2_ow1); % SE

m1_lc1_ow2  = mean(nbadcomps_exp1_lc1_2Hz_hc2_ow2,1);
m2_lc1_ow2  = mean(nbadcomps_exp2_lc1_2Hz_hc2_ow2,1);
se1_lc1_ow2 =  std(nbadcomps_exp1_lc1_2Hz_hc2_ow2); % SE
se2_lc1_ow2 =  std(nbadcomps_exp2_lc1_2Hz_hc2_ow2); % SE


m1_lc2_ow1  = mean(nbadcomps_exp1_lc2_2Hz_hc2_ow1,1);
m2_lc2_ow1  = mean(nbadcomps_exp2_lc2_2Hz_hc2_ow1,1);
se1_hc2_ow1 =  std(nbadcomps_exp1_lc2_2Hz_hc2_ow1); % alternative stderr()
se2_hc2_ow1 =  std(nbadcomps_exp2_lc2_2Hz_hc2_ow1); % SE

m1_hc2_ow2  = mean(nbadcomps_exp1_lc2_2Hz_hc2_ow2,1);
m2_hc2_ow2  = mean(nbadcomps_exp2_lc2_2Hz_hc2_ow2,1);
se1_hc2_ow2 =  std(nbadcomps_exp1_lc2_2Hz_hc2_ow2); % SE
se2_hc2_ow2 =  std(nbadcomps_exp2_lc2_2Hz_hc2_ow2); % SE


%% plot number of rejected ICs per threshold
subplot(6,2,3); hold on;
title('Scenes')
% individual subjects
for s = SUBJECTS
    x1 = nbadcomps_exp1_lc1_2Hz_hc2_ow1(s,:);
    hs1(s) = plot(THRESH2PLOT,x1,'marker','.','markersize',10);
    hs1(s).Color = [0.7,0.7,0.9, 0.9]; % R,G,B,Alpha
    hs1(s).MarkerFaceColor = [0.7,0.7,0.9];
    hs1(s).LineWidth = 0.7;
    if s == 12
        dummy = plot(THRESH2PLOT,x1,'marker','.','markersize',10);
        dummy.Color = [0.7,0.7,0.9, 0.9]; % R,G,B,Alpha
        dummy.MarkerFaceColor = [0.7,0.7,0.9];
        dummy.LineWidth = 0.7;
    end
end
h1 = plot(THRESH2PLOT,m1_lc1_ow1,'bo-','linewidth',1.2);
plot([THRESH2PLOT;THRESH2PLOT],[m1_lc1_ow1-se1_lc1_ow1;m1_lc1_ow1+se1_lc1_ow1],'b','linewidth',1.2);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
xlim([1.5 11.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:') % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
set(gca,'YTick',[0 15 30 45]);
lx = legend([h1 dummy],{'Average','Individ. subjects'})
set(lx,'box','off','Location','NorthEast')


subplot(6,2,4); hold on;
title('Reading')
% individual subjects
for s = SUBJECTS
    x2 = nbadcomps_exp2_lc1_2Hz_hc2_ow1(s,:);
    hs2(s) = plot(THRESH2PLOT,x2,'marker','.','markersize',10);
    hs2(s).Color = [0.9,0.7,0.7, 0.9]; % R,G,B,Alpha
    hs2(s).MarkerFaceColor = [0.9,0.7,0.7];
    hs2(s).LineWidth = 0.7;
    if s == 12
        dummy2 = plot(THRESH2PLOT,x1,'marker','.','markersize',10);
        dummy2.Color = [0.9,0.7,0.7, 0.9]; % R,G,B,Alpha
        dummy2.MarkerFaceColor = [0.9,0.7,0.7];
        dummy2.LineWidth = 0.7;
    end
end
h2 = plot(THRESH2PLOT,m2_lc1_ow1,'ro-','linewidth',1.2);

plot([THRESH2PLOT;THRESH2PLOT],[m2_lc1_ow1-se2_lc1_ow1;m2_lc1_ow1+se2_lc1_ow1],'r','linewidth',1.2);
xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
xlim([1.5 11.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:') % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs')
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
xlim([1.5 11.5])
ylim([0 NCHANS_EEG-1])
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
set(gca,'YTick',[0 15 30 45]);
l2 = legend([h2 dummy2],{'Average','Individ. subjects'})
set(l2,'box','off','Location','NorthEast')

%% by overweighting
subplot(6,2,5); hold on;
title('With/without overweighting')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh
MS = 4;
h10 = plot(THRESH2PLOT,m1_lc1_ow1,'ko-','linewidth',1.2);
h12 = plot(THRESH2PLOT,m1_lc1_ow2,'bo--','linewidth',1.2,'markersize',MS);
plot([THRESH2PLOT;THRESH2PLOT],[m1_lc1_ow1-se1_lc1_ow1;m1_lc1_ow1+se1_lc1_ow1],'k-','linewidth',1.2);
plot([THRESH2PLOT;THRESH2PLOT],[m1_lc1_ow2-se1_lc1_ow2;m1_lc1_ow2+se1_lc1_ow2],'b-','linewidth',1.2,'markersize',MS);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
xlim([1.5 11.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
set(gca,'YTick',[0 15 30 45]);
legend([h10 h12],{'basic','overweighted'},'box','off','Location','NorthEast');

subplot(6,2,6); hold on;
title('With/without overweighting')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh
MS = 4;
h11 = plot(THRESH2PLOT,m2_lc1_ow1,'ko-','linewidth',1.2);
h13 = plot(THRESH2PLOT,m2_lc1_ow2,'ro--','linewidth',1.2,'markersize',MS);
plot([THRESH2PLOT;THRESH2PLOT],[m2_lc1_ow1-se2_lc1_ow1;m2_lc1_ow1+se2_lc1_ow1],'k-','linewidth',1.2);
plot([THRESH2PLOT;THRESH2PLOT],[m2_lc1_ow2-se2_lc1_ow2;m2_lc1_ow2+se2_lc1_ow2],'r-','linewidth',1.2,'markersize',MS);
xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
%legend([h10 h11 h12 h13],{'Scenes: basic','Reading: basic','Scenes: overweighted','Reading: overweighted'},'box','off','Location','East');
xlim([1.5 11.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
set(gca,'YTick',[0 15 30 45]);
legend([h11 h13],{'basic','overweighted'},'box','off','Location','NorthEast');

% subplot(6,2,7); hold on; title('High-pass Filter')
%
% %% by HP-filter setting
% nbadcomps_exp1_hc2_ow2_filt1  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,4,OW_LEVEL,THRESH2PLOT));
% nbadcomps_exp1_hc2_ow2_filt2  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,5,OW_LEVEL,THRESH2PLOT));
% nbadcomps_exp1_hc2_ow2_filt3  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,10,OW_LEVEL,THRESH2PLOT));
% nbadcomps_exp1_hc2_ow2_filt4  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,15,OW_LEVEL,THRESH2PLOT));
% nbadcomps_exp1_hc2_ow2_filt5  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,16,OW_LEVEL,THRESH2PLOT));
%
% % subplot(6,2,6); hold on;
% % title('Number of rejected ICs by filter setting')
% % plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh
% % MS = 4;
% hf1 = plot(THRESH2PLOT,mean(nbadcomps_exp1_hc2_ow2_filt1),'bo-','linewidth',0.8);
% hf2 = plot(THRESH2PLOT,mean(nbadcomps_exp1_hc2_ow2_filt2),'bo-','linewidth',0.8);
% hf3 = plot(THRESH2PLOT,mean(nbadcomps_exp1_hc2_ow2_filt3),'bo-','linewidth',0.8);
% hf4 = plot(THRESH2PLOT,mean(nbadcomps_exp1_hc2_ow2_filt4),'bo-','linewidth',0.8);
% hf5 = plot(THRESH2PLOT,mean(nbadcomps_exp1_hc2_ow2_filt5),'bo-','linewidth',0.8);
%
% xlabel('Variance ratio threshold used')
% ylabel('N rejected ICs');
% %legend([h10 h11 h12 h13],{'Scenes: basic','Reading: basic','Scenes: overweighted','Reading: overweighted'},'box','off','Location','East');
% xlim([1.5 11.5])
% ylim([0 NCHANS_EEG-1])
% plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% set(gca,'XTick',t2plot_numbs);
% set(gca,'XTickLabel',t2plot_names);
% ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
% set(gca,'YTick',[0 15 30 45]);


%% CR, SP, Overcorrection by Threshold

% CR
metricsac_exp1_ow1_2Hz = squeeze(metric_sac(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricsac_exp2_ow1_2Hz = squeeze(metric_sac(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricsac_exp1_ow2_2Hz = squeeze(metric_sac(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));
metricsac_exp2_ow2_2Hz = squeeze(metric_sac(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));

% SP
metricsp_exp1_ow1_2Hz = squeeze(metric_sp(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricsp_exp2_ow1_2Hz = squeeze(metric_sp(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricsp_exp1_ow2_2Hz = squeeze(metric_sp(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));
metricsp_exp2_ow2_2Hz = squeeze(metric_sp(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));

% distortion
metricstim_exp1_ow1_2Hz = squeeze(metric_stim(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricstim_exp2_ow1_2Hz = squeeze(metric_stim(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,1,THRESH2PLOT));
metricstim_exp1_ow2_2Hz = squeeze(metric_stim(1,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));
metricstim_exp2_ow2_2Hz = squeeze(metric_stim(2,SUBJECTS,HC_LEVEL,LC_TO_PLOT,2,THRESH2PLOT));

% scenes
m1_cr_0  = mean(metricsac_exp1_ow1_2Hz,1);
m1_cr_1  = mean(metricsac_exp1_ow2_2Hz,1);
sd1_cr_0 =  std(metricsac_exp1_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd1_cr_1 =  std(metricsac_exp1_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));
m1_sp_0  = mean(metricsp_exp1_ow1_2Hz,1)
m1_sp_1  = mean(metricsp_exp1_ow2_2Hz,1)
sd1_sp_0 =  std(metricsp_exp1_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd1_sp_1 =  std(metricsp_exp1_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));
m1_di_0  = mean(metricstim_exp1_ow1_2Hz,1);
m1_di_1  = mean(metricstim_exp1_ow2_2Hz,1);
sd1_di_0 =  std(metricstim_exp1_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd1_di_1 =  std(metricstim_exp1_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));

% reading
m2_cr_0  = mean(metricsac_exp2_ow1_2Hz,1);
m2_cr_1  = mean(metricsac_exp2_ow2_2Hz,1);
sd2_cr_0 =  std(metricsac_exp2_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd2_cr_1 =  std(metricsac_exp2_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));
m2_sp_0  = mean(metricsp_exp2_ow1_2Hz,1)
m2_sp_1  = mean(metricsp_exp2_ow2_2Hz,1)
sd2_sp_0 =  std(metricsp_exp2_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd2_sp_1 =  std(metricsp_exp2_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));
m2_di_0  = mean(metricstim_exp2_ow1_2Hz,1);
m2_di_1  = mean(metricstim_exp2_ow2_2Hz,1);
sd2_di_0 =  std(metricstim_exp2_ow1_2Hz,0,1) / sqrt(length(SUBJECTS));
sd2_di_1 =  std(metricstim_exp2_ow2_2Hz,0,1) / sqrt(length(SUBJECTS));

% f12 = figure('name','Effect of rejection threshold (at 2 Hz)')

% % scenes
subplot(6,2,9); hold on; title('Corneoretinal artifact');
% h1 = plot([THRESH2PLOT],m1_cr_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT],m1_cr_1,'bo-','linewidth',LW)
% h1 = plot([THRESH2PLOT], mean(metricsac_exp1_ow1_2Hz,1),'bx:','linewidth',LW)
% h2 = plot([THRESH2PLOT], mean(metricsac_exp1_ow2_2Hz,1),'bo-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m1_cr_0-sd1_cr_0;m1_cr_0+sd1_cr_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m1_cr_1-sd1_cr_1;m1_cr_1+sd1_cr_1],'b-','linewidth',LW2);
% lxx = legend([h1 h2],{'basic','overweighted'},'box','off')
set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([1.5 11.5])
ylim([-2 4])
xlabel('Threshold')
ylabel('Frontal lateralization [µV]')
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(6,2,11); hold on; title('Spike potential');
% h1 = plot([THRESH2PLOT],m1_sp_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT],m1_sp_1,'bo-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m1_sp_0-sd1_sp_0;m1_sp_0+sd1_sp_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m1_sp_1-sd1_sp_1;m1_sp_1+sd1_sp_1],'b-','linewidth',LW2);
set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([1.5 11.5])
ylim([0 2.1])
xlabel('Threshold')
ylabel('GFP [µV]')

ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(6,2,[7]); hold on; title('ERP distortion');
h2 = plot([THRESH2PLOT],m1_di_1,'bo-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m1_di_0-sd1_di_0;m1_di_0+sd1_di_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m1_di_1-sd1_di_1;m1_di_1+sd1_di_1],'b-','linewidth',LW2);set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
xlim([1.5 11.5])
ylim([0.35 1.75])
xlabel('Variance ratio threshold used')
ylabel('GFP of ERP distortion [µV]')

ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

% reading
subplot(6,2,10); hold on; title('Corneoretinal artifact');
% h1 = plot([THRESH2PLOT],m2_cr_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT],m2_cr_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m2_cr_0-sd2_cr_0;m2_cr_0+sd2_cr_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m2_cr_1-sd2_cr_1;m2_cr_1+sd2_cr_1],'r-','linewidth',LW2);
set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([1.5 11.5])
ylim([-2 4])
xlabel('Threshold')
ylabel('Frontal lateralization [µV]')

ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(6,2,12); hold on; title('Spike potential');
% h1 = plot([THRESH2PLOT],m2_sp_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT],m2_sp_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m2_sp_0-sd2_sp_0;m2_sp_0+sd2_sp_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m2_sp_1-sd2_sp_1;m2_sp_1+sd2_sp_1],'r-','linewidth',LW2);
set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([1.5 11.5])
ylim([0 2.1])
xlabel('Threshold')
ylabel('GFP [µV]')

ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(6,2,[8]); hold on; title('ERP distortion');
% h1 = plot([THRESH2PLOT],m2_di_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT],m2_di_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT;THRESH2PLOT],[m2_di_0-sd2_di_0;m2_di_0+sd2_di_0],'k-','linewidth',LW2);
plot([THRESH2PLOT;THRESH2PLOT],[m2_di_1-sd2_di_1;m2_di_1+sd2_di_1],'r-','linewidth',LW2);set(gca,'Xtick',THRESH2PLOT);
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);
xlim([1.5 11.5])
ylim([0.35 1.75])
xlabel('Variance ratio threshold used')
ylabel('GFP of ERP distortion [µV]')
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

% export figure
if EXPORTFIG
    set(f_thresh,'Position',SCREENSIZE)
    export_fig(f_thresh,'figs_2019_R1/fig_thresholds_v9_ProofCorrection_scale0_25','-eps','-pdf','-transparent','-painters')
end