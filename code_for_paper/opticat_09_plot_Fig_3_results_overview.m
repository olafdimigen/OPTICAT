%% PLOT RESULTS OF "OPTICA" STUDY
% olaf.dimigen@hu-berlin.de, 2018
clear, rng('shuffle'), close all

addpath('M:/Dropbox/_subfunc_master') % functions: stderr und ci ...
addpath M:/Dropbox/eeglab14_1_1b; eeglab; close
addpath M:/Dropbox/_subfunc_master/export_fig_2018 % export_fig
addpath M:/Dropbox/_subfunc_master/ufo/newplot % putnewaxes by Ulrich R.

LC  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30];
HC  = [40 100];
THR = [0.5:0.1:1.6]; % thresholds, 1.6 = outlier-detection

SUBJECTS   = 1:12;
EXPORTFIG  = 1; % export figures as .eps?
SCREENSIZE = [10 10 1000 1000]; % size for figure export (pixel)

%% load input data

% with GFP zero point
load Z:/OPTICA/results/metrics.mat 
load Z:/OPTICA/results/opticaresults.mat opticaresults LOWCUTOFFS HICUTOFFS THRESHOLDS
load Z:/OPTICA/results/results_for_plotting.mat SCENE* READ*
load Z:/OPTICA/results/chanlocs.mat chanlocs_reading chanlocs_scenes

%% -- SUMMARZ: Aggregated results at threshold of 1.1 (filters x overweighting)
THRESHOLD2PLOT = 7 % = threshold of 1.1
FILTERLABELS   = {'No','.1','.25','.5','.75','1','1.5','2','2.5','3','3.5','4','5','7.5','10','12.5','15','20','25','30','MSEC'} % passband edges
LW             = 1.0; % lines
LW2            = 0.5; % error bars
MS             =   5; % marker size
FSA            =  12; % font size axis

% plot colors (lines, markers) etc.
COLOR_1    = [204,76,2]    ./255; % 
COLOR_3    = [254,153,41]  ./255  % 
COLOR_2    = [34,94,168]   ./255  %
COLOR_4    = [65,182,196]  ./255  %
COLOR_MSEC = [120,198,121] ./255; % 

LINESTYLE_1 = '-';
LINESTYLE_3 = ':';
LINESTYLE_2 = '-';
LINESTYLE_4 = ':';

MARKER_1    = '^';
MARKER_3    = '^';
MARKER_2    = 'o';
MARKER_4    = 'o';
MARKER_MSEC = 'o';

FACEALPHA     = 0.1
ANGLE_XLABELS = 60; % angle of x axis labels

% CR: 40 Hz
metricsac_exp1_hc1_ow1 = squeeze(metric_sac(1,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricsac_exp2_hc1_ow1 = squeeze(metric_sac(2,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricsac_exp1_hc1_ow2 = squeeze(metric_sac(1,SUBJECTS,1,:,2,THRESHOLD2PLOT));
metricsac_exp2_hc1_ow2 = squeeze(metric_sac(2,SUBJECTS,1,:,2,THRESHOLD2PLOT));
% CR: 100 Hz
metricsac_exp1_hc2_ow1 = squeeze(metric_sac(1,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricsac_exp2_hc2_ow1 = squeeze(metric_sac(2,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricsac_exp1_hc2_ow2 = squeeze(metric_sac(1,SUBJECTS,2,:,2,THRESHOLD2PLOT));
metricsac_exp2_hc2_ow2 = squeeze(metric_sac(2,SUBJECTS,2,:,2,THRESHOLD2PLOT));

% SP: 40 Hz
metricsp_exp1_hc1_ow1 = squeeze(metric_sp(1,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricsp_exp2_hc1_ow1 = squeeze(metric_sp(2,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricsp_exp1_hc1_ow2 = squeeze(metric_sp(1,SUBJECTS,1,:,2,THRESHOLD2PLOT));
metricsp_exp2_hc1_ow2 = squeeze(metric_sp(2,SUBJECTS,1,:,2,THRESHOLD2PLOT));
% SP: 100 Hz
metricsp_exp1_hc2_ow1 = squeeze(metric_sp(1,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricsp_exp2_hc2_ow1 = squeeze(metric_sp(2,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricsp_exp1_hc2_ow2 = squeeze(metric_sp(1,SUBJECTS,2,:,2,THRESHOLD2PLOT));
metricsp_exp2_hc2_ow2 = squeeze(metric_sp(2,SUBJECTS,2,:,2,THRESHOLD2PLOT));
 
% new: add baseline GFP for comparison

% SP: 40 Hz
metricbasegfp_exp1_hc1_ow1 = squeeze(metric_baseGFP(1,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricbasegfp_exp2_hc1_ow1 = squeeze(metric_baseGFP(2,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricbasegfp_exp1_hc1_ow2 = squeeze(metric_baseGFP(1,SUBJECTS,1,:,2,THRESHOLD2PLOT));
metricbasegfp_exp2_hc1_ow2 = squeeze(metric_baseGFP(2,SUBJECTS,1,:,2,THRESHOLD2PLOT));
% SP: 100 Hz
metricbasegfp_exp1_hc2_ow1 = squeeze(metric_baseGFP(1,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricbasegfp_exp2_hc2_ow1 = squeeze(metric_baseGFP(2,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricbasegfp_exp1_hc2_ow2 = squeeze(metric_baseGFP(1,SUBJECTS,2,:,2,THRESHOLD2PLOT));
metricbasegfp_exp2_hc2_ow2 = squeeze(metric_baseGFP(2,SUBJECTS,2,:,2,THRESHOLD2PLOT));

% Overcorrection: 40 Hz
metricstim_exp1_hc1_ow1 = squeeze(metric_stim(1,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricstim_exp2_hc1_ow1 = squeeze(metric_stim(2,SUBJECTS,1,:,1,THRESHOLD2PLOT));
metricstim_exp1_hc1_ow2 = squeeze(metric_stim(1,SUBJECTS,1,:,2,THRESHOLD2PLOT));
metricstim_exp2_hc1_ow2 = squeeze(metric_stim(2,SUBJECTS,1,:,2,THRESHOLD2PLOT));
% Overcorrection: 100 Hz
metricstim_exp1_hc2_ow1 = squeeze(metric_stim(1,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricstim_exp2_hc2_ow1 = squeeze(metric_stim(2,SUBJECTS,2,:,1,THRESHOLD2PLOT));
metricstim_exp1_hc2_ow2 = squeeze(metric_stim(1,SUBJECTS,2,:,2,THRESHOLD2PLOT));
metricstim_exp2_hc2_ow2 = squeeze(metric_stim(2,SUBJECTS,2,:,2,THRESHOLD2PLOT));

%% Aggregate: Scenes: CR
m1_cr_hc1_ow1  = mean(metricsac_exp1_hc1_ow1,1);
m1_cr_hc1_ow2  = mean(metricsac_exp1_hc1_ow2,1);
sd1_cr_hc1_ow1 =  std(metricsac_exp1_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_cr_hc1_ow2 =  std(metricsac_exp1_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m1_cr_hc2_ow1  = mean(metricsac_exp1_hc2_ow1,1);
m1_cr_hc2_ow2  = mean(metricsac_exp1_hc2_ow2,1);
sd1_cr_hc2_ow1 =  std(metricsac_exp1_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_cr_hc2_ow2 =  std(metricsac_exp1_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

%% Aggregate: Scenes: SP
m1_sp_hc1_ow1  = mean(metricsp_exp1_hc1_ow1,1);
m1_sp_hc1_ow2  = mean(metricsp_exp1_hc1_ow2,1);
sd1_sp_hc1_ow1 =  std(metricsp_exp1_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_sp_hc1_ow2 =  std(metricsp_exp1_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m1_sp_hc2_ow1  = mean(metricsp_exp1_hc2_ow1,1);
m1_sp_hc2_ow2  = mean(metricsp_exp1_hc2_ow2,1);
sd1_sp_hc2_ow1 =  std(metricsp_exp1_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_sp_hc2_ow2 =  std(metricsp_exp1_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

%% Aggregate: Scenes: Base GFP
% SP: 40 Hz
m1_basegfp_hc1_ow1  = mean(metricbasegfp_exp1_hc1_ow1,1);
m1_basegfp_hc1_ow2  = mean(metricbasegfp_exp1_hc1_ow2,1);
sd1_basegfp_hc1_ow1 = std(metricbasegfp_exp1_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_basegfp_hc1_ow2 = std(metricbasegfp_exp1_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

% SP: 100 Hz
m1_basegfp_hc2_ow1  = mean(metricbasegfp_exp1_hc2_ow1,1);
m1_basegfp_hc2_ow2  = mean(metricbasegfp_exp1_hc2_ow2,1);
sd1_basegfp_hc2_ow1 = std(metricbasegfp_exp1_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_basegfp_hc2_ow2 = std(metricbasegfp_exp1_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

% note: baseline GFP is extremely similar for different levels of the high-pass filter (etc)
% for clarity, just plot grand mean of the baseline GFP: 
m1_basegfp_grandavg = (mean(m1_basegfp_hc1_ow1) + mean(m1_basegfp_hc1_ow2) + mean(m1_basegfp_hc2_ow1) + mean(m1_basegfp_hc2_ow2)) ./4;

%% Aggregate: Scenes: Distortion
m1_di_hc1_ow1  = mean(metricstim_exp1_hc1_ow1,1);
m1_di_hc1_ow2  = mean(metricstim_exp1_hc1_ow2,1);
sd1_di_hc1_ow1 =  std(metricstim_exp1_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_di_hc1_ow2 =  std(metricstim_exp1_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m1_di_hc2_ow1  = mean(metricstim_exp1_hc2_ow1,1);
m1_di_hc2_ow2  = mean(metricstim_exp1_hc2_ow2,1);
sd1_di_hc2_ow1 =  std(metricstim_exp1_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd1_di_hc2_ow2 =  std(metricstim_exp1_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

%% Aggregate: Reading: CR
m2_cr_hc1_ow1  = mean(metricsac_exp2_hc1_ow1,1);
m2_cr_hc1_ow2  = mean(metricsac_exp2_hc1_ow2,1);
sd2_cr_hc1_ow1 =  std(metricsac_exp2_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_cr_hc1_ow2 =  std(metricsac_exp2_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m2_cr_hc2_ow1  = mean(metricsac_exp2_hc2_ow1,1);
m2_cr_hc2_ow2  = mean(metricsac_exp2_hc2_ow2,1);
sd2_cr_hc2_ow1 =  std(metricsac_exp2_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_cr_hc2_ow2 =  std(metricsac_exp2_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

%% Aggregate: Reading: SP
m2_sp_hc1_ow1  = mean(metricsp_exp2_hc1_ow1,1);
m2_sp_hc1_ow2  = mean(metricsp_exp2_hc1_ow2,1);
sd2_sp_hc1_ow1 =  std(metricsp_exp2_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_sp_hc1_ow2 =  std(metricsp_exp2_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m2_sp_hc2_ow1  = mean(metricsp_exp2_hc2_ow1,1);
m2_sp_hc2_ow2  = mean(metricsp_exp2_hc2_ow2,1);
sd2_sp_hc2_ow1 =  std(metricsp_exp2_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_sp_hc2_ow2 =  std(metricsp_exp2_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

%% Aggregate: Reading: Base GFP
% SP: 40 Hz
m2_basegfp_hc1_ow1  = mean(metricbasegfp_exp2_hc1_ow1,1);
m2_basegfp_hc1_ow2  = mean(metricbasegfp_exp2_hc1_ow2,1);
sd2_basegfp_hc1_ow1 = std(metricbasegfp_exp2_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_basegfp_hc1_ow2 = std(metricbasegfp_exp2_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

% SP: 100 Hz
m2_basegfp_hc2_ow1  = mean(metricbasegfp_exp2_hc2_ow1,1);
m2_basegfp_hc2_ow2  = mean(metricbasegfp_exp2_hc2_ow2,1);
sd2_basegfp_hc2_ow1 = std(metricbasegfp_exp2_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_basegfp_hc2_ow2 = std(metricbasegfp_exp2_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

m2_basegfp_grandavg = (mean(m2_basegfp_hc1_ow1) + mean(m2_basegfp_hc1_ow2) + mean(m2_basegfp_hc2_ow1) + mean(m2_basegfp_hc2_ow2)) ./4;

%% Aggregate: Reading: Distortion
m2_di_hc1_ow1  = mean(metricstim_exp2_hc1_ow1,1);
m2_di_hc1_ow2  = mean(metricstim_exp2_hc1_ow2,1);
sd2_di_hc1_ow1 =  std(metricstim_exp2_hc1_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_di_hc1_ow2 =  std(metricstim_exp2_hc1_ow2,0,1) / sqrt(length(SUBJECTS));

m2_di_hc2_ow1  = mean(metricstim_exp2_hc2_ow1,1);
m2_di_hc2_ow2  = mean(metricstim_exp2_hc2_ow2,1);
sd2_di_hc2_ow1 =  std(metricstim_exp2_hc2_ow1,0,1) / sqrt(length(SUBJECTS));
sd2_di_hc2_ow2 =  std(metricstim_exp2_hc2_ow2,0,1) / sqrt(length(SUBJECTS));

% BESA: scenes
m1_cr_msec  = mean(metric_sac_msec(1,:),2);
sd1_cr_msec =  std(metric_sac_msec(1,:),0,2) / sqrt(length(SUBJECTS));
m1_sp_msec  = mean(metric_sp_msec(1,:),2);
sd1_sp_msec =  std(metric_sp_msec(1,:),0,2) / sqrt(length(SUBJECTS));
m1_di_msec  = mean(metric_stim_msec(1,:),2);
sd1_di_msec =  std(metric_stim_msec(1,:),0,2) / sqrt(length(SUBJECTS));

% BESA: reading
m2_cr_msec  = mean(metric_sac_msec(2,:),2);
sd2_cr_msec =  std(metric_sac_msec(2,:),0,2) / sqrt(length(SUBJECTS));
m2_sp_msec  = mean(metric_sp_msec(2,:),2);
sd2_sp_msec =  std(metric_sp_msec(2,:),0,2) / sqrt(length(SUBJECTS));
m2_di_msec  = mean(metric_stim_msec(2,:),2);
sd2_di_msec =  std(metric_stim_msec(2,:),0,2) / sqrt(length(SUBJECTS));


%% FIGURE
f1 = figure

%% scenes: CR
subplot(3,2,1); hold on; title('Corneoretinal');
plot([1 22],[0 0],'k:') % voltage zero

% % 40 Hz: Lines: error bars
plot([1:20;1:20],[m1_cr_hc1_ow1-sd1_cr_hc1_ow1;m1_cr_hc1_ow1+sd1_cr_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m1_cr_hc1_ow2-sd1_cr_hc1_ow2;m1_cr_hc1_ow2+sd1_cr_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m1_cr_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle','-','color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m1_cr_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

%100 Hz: Lines: error bars
plot([1:20;1:20],[m1_cr_hc2_ow1-sd1_cr_hc2_ow1;m1_cr_hc2_ow1+sd1_cr_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m1_cr_hc2_ow2-sd1_cr_hc2_ow2;m1_cr_hc2_ow2+sd1_cr_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m1_cr_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2)
h4 = plot(1:20,m1_cr_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% MSEC
plot([22 22],[m1_cr_msec-sd1_cr_msec;m1_cr_msec+sd1_cr_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m1_cr_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point

l = legend([h1 h3 h2 h4 mm],{'HP 40 Hz | basic','HP 40 Hz | overweighted','HP 100 Hz | basic','HP 100 Hz | overweighted','MSEC'},'box','off')

set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([-7 27])
ylabel('CRA: Frontal lateraliz. [µV]')
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
set(gca,'fontsize',FSA)
% xtickangle(ANGLE_XLABELS)
% make plot prettier
[x_pna1, y_pna1] = putnewaxes([1 22],[-5 25],[-8],[0])
% xtickangle(x_pna1,ANGLE_XLABELS)
set(x_pna1,'Xtick',[1:20 22]);
set(x_pna1,'Xticklabels',FILTERLABELS);

%% reading: CR
subplot(3,2,2); hold on; title('Corneoretinal')
plot([1 22],[0 0],'k:') % voltage zero

% 40 Hz: Lines: error bars
plot([1:20;1:20],[m2_cr_hc1_ow1-sd2_cr_hc1_ow1;m2_cr_hc1_ow1+sd2_cr_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m2_cr_hc1_ow2-sd2_cr_hc1_ow2;m2_cr_hc1_ow2+sd2_cr_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m2_cr_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle',LINESTYLE_1,'color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m2_cr_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

% 100 Hz: Lines: error bars
plot([1:20;1:20],[m2_cr_hc2_ow1-sd2_cr_hc2_ow1;m2_cr_hc2_ow1+sd2_cr_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m2_cr_hc2_ow2-sd2_cr_hc2_ow2;m2_cr_hc2_ow2+sd2_cr_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m2_cr_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2)
h4 = plot(1:20,m2_cr_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% MSEC
plot([22 22],[m2_cr_msec-sd2_cr_msec;m2_cr_msec+sd2_cr_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m2_cr_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point

%l = legend([h1 h3 h2 h4 mm],{'HP 40 Hz | basic','HP 40 Hz | overweighted','HP 100 Hz | basic','HP 100 Hz | overweighted','MSEC'},'box','off')

set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([-7 27])
ylabel('CRA: Frontal lateraliz. [µV]')
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
set(gca,'fontsize',FSA)
% xtickangle(ANGLE_XLABELS)
% make plot prettier
[x_pna2, y_pna2] = putnewaxes([1 22],[-5 25],[-8],[0])
% xtickangle(x_pna2,ANGLE_XLABELS)
set(x_pna2,'Xtick',[1:20 22]);
set(x_pna2,'Xticklabels',FILTERLABELS);


%% Scenes: SP
subplot(3,2,3); hold on; title('Spike potential')

% 40 Hz: Lines: error bars
plot([1:20;1:20],[m1_sp_hc1_ow1-sd1_sp_hc1_ow1;m1_sp_hc1_ow1+sd1_sp_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m1_sp_hc1_ow2-sd1_sp_hc1_ow2;m1_sp_hc1_ow2+sd1_sp_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m1_sp_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle',LINESTYLE_1,'color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m1_sp_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

% baseline GFP, simplified
h100 = plot([1 22],[m1_basegfp_grandavg m1_basegfp_grandavg],'k:')

% 100 Hz: Lines: error bars
plot([1:20;1:20],[m1_sp_hc2_ow1-sd1_sp_hc2_ow1;m1_sp_hc2_ow1+sd1_sp_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m1_sp_hc2_ow2-sd1_sp_hc2_ow2;m1_sp_hc2_ow2+sd1_sp_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m1_sp_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2)
h4 = plot(1:20,m1_sp_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% baseline GFP, simplified
h100 = plot([1 22],[m1_basegfp_grandavg m1_basegfp_grandavg],'k:')

% MSEC
plot([22 22],[m1_sp_msec-sd1_sp_msec;m1_sp_msec+sd1_sp_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m1_sp_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point

set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([0 6.3])
ylabel('SP amplitude, GFP [µV]')
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
% xtickangle(ANGLE_XLABELS)
set(gca,'fontsize',FSA)
% xtickangle(ANGLE_XLABELS)
% make plot prettier
[x_pna3, y_pna3] = putnewaxes([1 22],[0 6],[-1],[0])
% xtickangle(x_pna3,ANGLE_XLABELS)
set(x_pna3,'Xtick',[1:20 22]);
set(x_pna3,'Xticklabels',FILTERLABELS);

%% reading: SP
subplot(3,2,4); hold on; title('Spike potential')

% 40 Hz: Lines: error bars
plot([1:20;1:20],[m2_sp_hc1_ow1-sd2_sp_hc1_ow1;m2_sp_hc1_ow1+sd2_sp_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m2_sp_hc1_ow2-sd2_sp_hc1_ow2;m2_sp_hc1_ow2+sd2_sp_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m2_sp_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle',LINESTYLE_1,'color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m2_sp_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

% baseline GFP, simplified
h200 = plot([1 22],[m2_basegfp_grandavg m2_basegfp_grandavg],'k:')

% 100 Hz: Lines: error bars
plot([1:20;1:20],[m2_sp_hc2_ow1-sd2_sp_hc2_ow1;m2_sp_hc2_ow1+sd2_sp_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m2_sp_hc2_ow2-sd2_sp_hc2_ow2;m2_sp_hc2_ow2+sd2_sp_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m2_sp_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2)
h4 = plot(1:20,m2_sp_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% baseline GFP, simplified
h200 = plot([1 22],[m2_basegfp_grandavg m2_basegfp_grandavg],'k:')

% MSEC
plot([22 22],[m2_sp_msec-sd2_sp_msec;m2_sp_msec+sd2_sp_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m2_sp_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point

set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([0 6.3])
ylabel('SP amplitude, GFP [µV]')
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
% xtickangle(ANGLE_XLABELS)
set(gca,'fontsize',FSA)
% make plot prettier
[x_pna4, y_pna4] = putnewaxes([1 22],[0 6],[-1],[0])
% xtickangle(x_pna4,ANGLE_XLABELS)
set(x_pna4,'Xtick',[1:20 22]);
set(x_pna4,'Xticklabels',FILTERLABELS);



%% scenes: distortion
subplot(3,2,5); hold on; title('Overcorrection')

% 40 Hz: Lines: error bars
plot([1:20;1:20],[m1_di_hc1_ow1-sd1_di_hc1_ow1;m1_di_hc1_ow1+sd1_di_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m1_di_hc1_ow2-sd1_di_hc1_ow2;m1_di_hc1_ow2+sd1_di_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m1_di_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle',LINESTYLE_1,'color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m1_di_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

% 100 Hz: Lines: error bars
plot([1:20;1:20],[m1_di_hc2_ow1-sd1_di_hc2_ow1;m1_di_hc2_ow1+sd1_di_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m1_di_hc2_ow2-sd1_di_hc2_ow2;m1_di_hc2_ow2+sd1_di_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m1_di_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2);
h4 = plot(1:20,m1_di_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% MSEC
plot([22 22],[m1_di_msec-sd1_di_msec;m1_di_msec+sd1_di_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m1_di_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point

set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([0 3.8])
xlabel('High-pass filter (Hz)')
ylabel('ERP distortion, GFP [µV]')
% axis square, box on
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
% xtickangle(ANGLE_XLABELS)
set(gca,'fontsize',FSA)

% make plot more pretty
[x_pna5, y_pna5] = putnewaxes([1 22],[0 3.5],[-1],[0])
% % xtickangle(x_pna5,90)
set(x_pna5,'Xtick',[1:20 22]);
set(x_pna5,'Xticklabels',FILTERLABELS);


%% distortion: reading
subplot(3,2,6); hold on; title('Overcorrection')

% 40 Hz: Lines: error bars
plot([1:20;1:20],[m2_di_hc1_ow1-sd2_di_hc1_ow1;m2_di_hc1_ow1+sd2_di_hc1_ow1],'linestyle','-','color',COLOR_1,'linewidth',LW2);
plot([1:20;1:20],[m2_di_hc1_ow2-sd2_di_hc1_ow2;m2_di_hc1_ow2+sd2_di_hc1_ow2],'linestyle','-','color',COLOR_3,'linewidth',LW2);
% 40 Hz: Lines
h1 = plot(1:20,m2_di_hc1_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_1,'linestyle',LINESTYLE_1,'color',COLOR_1,'marker',MARKER_1)
h3 = plot(1:20,m2_di_hc1_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_3,'linestyle',LINESTYLE_3,'color',COLOR_3,'marker',MARKER_3)

% 100 Hz: Lines: error bars
plot([1:20;1:20],[m2_di_hc2_ow1-sd2_di_hc2_ow1;m2_di_hc2_ow1+sd2_di_hc2_ow1],'linestyle','-','color',COLOR_2,'linewidth',LW2);
plot([1:20;1:20],[m2_di_hc2_ow2-sd2_di_hc2_ow2;m2_di_hc2_ow2+sd2_di_hc2_ow2],'linestyle','-','color',COLOR_4,'linewidth',LW2);
% 100 Hz: Lines
h2 = plot(1:20,m2_di_hc2_ow1,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_2,'linestyle',LINESTYLE_2,'color',COLOR_2,'marker',MARKER_2)
h4 = plot(1:20,m2_di_hc2_ow2,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_4,'linestyle',LINESTYLE_4,'color',COLOR_4,'marker',MARKER_4)

% MSEC
plot([22 22],[m2_di_msec-sd2_di_msec;m2_di_msec+sd2_di_msec],'color',COLOR_MSEC,'linewidth',LW2); % SD
mm = plot([22],[m2_di_msec],'color',COLOR_MSEC,'linewidth',LW,'markersize',MS,'MarkerFaceColor',COLOR_MSEC,'marker',MARKER_MSEC); % data point


set(gca,'Xtick',[1:20 22] );
set(gca,'Xticklabels',FILTERLABELS);
xlim([0.2 22.8])
ylim([0 2.0])
xlabel('High-pass filter (Hz)')
ylabel('ERP distortion, GFP [µV]')
ylimits = ylim; % plot([9.5 9.5],[ylimits(1) ylimits(2)],'k:')
% xtickangle(ANGLE_XLABELS)
set(gca,'fontsize',FSA)
% make plot prettier
[x_pna6, y_pna6] = putnewaxes([1 22],[0 2.0],[-1],[0])
% xtickangle(x_pna6,ANGLE_XLABELS)
set(x_pna6,'Xtick',[1:20 22]);
set(x_pna6,'Xticklabels',FILTERLABELS);

%% title for whole figure
% [axt,ht]=suplabel('Scenes                                      Reading','t');
% set(ht,'FontSize',18)
% [axty,hty]=suplabel('READING     /     SCENES','y');
% set(hty,'FontSize',18)

set(f1,'Position',SCREENSIZE)

% % export figure
if EXPORTFIG
    export_fig(f1,'figs_2018/figure_metric_bigsummary_v3_R2_withBaselineGFP','-eps','-transparent','-painters')
end


% #################################################################################################################

%% FIGURE: distribution of variance ratios
all_vr1 = []; all_vr2 = []; all_vr12 = [];
for exp = 1:2
    for s = 1:12
        
        vr = opticaresults(exp).subj(s).hc(2).lc(8).ow(1).varratio50; % 8 = 2 Hz
        % store
        if exp == 1
            all_vr1 = [all_vr1; vr];
        else
            all_vr2 = [all_vr2; vr];
        end
        all_vr12 = [all_vr12; vr];
    end
end

MAX = 1.5 % summarize all values equal to or above this
all_vr1(all_vr1 > MAX)  = MAX;
all_vr2(all_vr2 > MAX)  = MAX;
all_vr12(all_vr12 > MAX) = MAX;

edges = [0:0.1:10];
n1 = histc(all_vr1  ,edges);
n2 = histc(all_vr2  ,edges);
n3 = histc(all_vr12 ,edges);

% divide by 45 (normalize)
n1 = n1/45;
n2 = n2/45;
n3 = n3/45;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X-Ticks and X-Tick labels
t2plot_numbs = [1 3 5 7 9 11 12];
t2plot_names = {'0.5','0.7','0.9','1.1','1.3','1.5','auto'};


%% plot figure
figure;
subplot(4,2,1); title('Scenes'); hold on
plot([1 1],[0 4.5],'r'); h1 = bar(edges,n1); xlim([0 MAX+0.05]);
ylabel('Number of ICs')
xlabel('Variance ratio threshold used')
set(gca,'XTick',[0.5:0.1:1.6]);
set(gca,'XTickLabel',{'0','1','2','3','>4'});
h1.EdgeColor = 'none'; h1.BarWidth = 0.75;
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names)
ylim([0 3])

subplot(4,2,2); title('Reading'); hold on
plot([1 1],[0 4.5],'r'); h2 = bar(edges,n2); xlim([0 MAX+0.05]);
ylabel('Number of ICs')
xlabel('Variance ratio threshold used')
% set(gca,'XTick',0:4);
% set(gca,'XTickLabel',{'0','1','2','3','>4'});
h2.EdgeColor = 'none'; h2.BarWidth = 0.75;
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names)
ylim([0 3])

% subplot(1,3,3); title('Both'); hold on
% plot([1 1],[0 4.5],'r'); h3 = bar(edges,n3); xlim([0 MAX+0.05]);
% ylabel('Number of ICs')
% xlabel('Saccade/fixation Variance Ratio')
% set(gca,'XTick',0:4);
% set(gca,'XTickLabel',{'0','1','2','3','>4'});
% h1.EdgeColor = 'none'; h1.BarWidth = 0.75;
% h2.EdgeColor = 'none'; h2.BarWidth = 0.75;
% h3.EdgeColor = 'none'; h3.BarWidth = 0.75;


% #################################################################################################################


%% FIGURE: N REJECTED IC (AS FUNCTION OF THRESHOLD)

% NOTE: size(nbadcomps) = exp x subjects x filter x version x thresholds
% lc level 8 = 2 Hz

nbadcomps_exp1_2Hz_hc2_ow1  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,8,1,:));
nbadcomps_exp1_2Hz_hc2_ow2  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,8,2,:));

nbadcomps_exp2_2Hz_hc2_ow1  = squeeze(nbadcomps(2,SUBJECTS,HC_LEVEL,8,1,:));
nbadcomps_exp2_2Hz_hc2_ow2  = squeeze(nbadcomps(2,SUBJECTS,HC_LEVEL,8,2,:));

m1_ow1  = mean(nbadcomps_exp1_2Hz_hc2_ow1,1);
m2_ow1  = mean(nbadcomps_exp2_2Hz_hc2_ow1,1);
se1_ow1 = stderr(nbadcomps_exp1_2Hz_hc2_ow1); % standard error
se2_ow1 = stderr(nbadcomps_exp2_2Hz_hc2_ow1); % standard error

m1_ow2  = mean(nbadcomps_exp1_2Hz_hc2_ow2,1);
m2_ow2  = mean(nbadcomps_exp2_2Hz_hc2_ow2,1);
se1_ow2 = stderr(nbadcomps_exp1_2Hz_hc2_ow2); % standard error
se2_ow2 = stderr(nbadcomps_exp2_2Hz_hc2_ow2); % standard error

NCHANS_EEG = 46;

subplot(4,2,3); hold on;
title('Scenes')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh

h1 = plot(1:length(THRESHOLDS),m1_ow1,'bo-','linewidth',1.2);
h2 = plot(1:length(THRESHOLDS),m2_ow1,'ro-','linewidth',1.2);
plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m1_ow1-se1_ow1;m1_ow1+se1_ow1],'b','linewidth',1.2);
plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m2_ow1-se2_ow1;m2_ow1+se2_ow1],'r','linewidth',1.2);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
legend([h1 h2],{'Scenes','Reading'},'box','off','Location','East');
xlim([0.5 12.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);


subplot(4,2,4); hold on;
title('Individuals')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh

for s = SUBJECTS
    
    x1 = nbadcomps_exp1_2Hz_hc2_ow1(s,:);
    x2 = nbadcomps_exp2_2Hz_hc2_ow1(s,:);
    hs1 = plot(1:length(THRESHOLDS),x1,'marker','.');
    hs2 = plot(1:length(THRESHOLDS),x2,'marker','.');
    
    % randomize line color, width, transparency
    r1 = rand./4;
    r2 = rand./4;
    hs1.Color = [r1*2,r2*2,0.75+r1,0.3+r1*2]; % R,G,B,Alpha
    hs1.MarkerFaceColor = [r1*2,r2*2,0.75+r1];% R,G,B
    hs1.LineWidth = 0.6+r1*4;
    hs2.Color = [0.75+r2,r1*2,r2*2,0.3+r2*2];
    hs2.MarkerFaceColor = [0.75+r2,r1*2,r2*2];
    hs2.LineWidth = 0.6+r2*4;
    
end

% if s == 1
legend([h1 h2],{'Scenes','Reading'},'box','off'); %'Location','East'
% end

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs')

% legend([h1 h2],{'Scenes','Reading'},'box','off');
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
xlim([0.5 12.5])
ylim([0 NCHANS_EEG-1])

% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);

% h3 = plot(1:length(THRESHOLDS),m1c,'bo:','linewidth',1.2); % owion 3 for comparison
% h4 = plot(1:length(THRESHOLDS),m2c,'ro:','linewidth',1.2);
% errorbar(m2,sd1,'bo-');
% errorbar(m2,sd2,'ro-');

subplot(4,2,5); hold on;
title('Number of rejected ICs: with/without extra SPs')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh
MS = 4;
h10 = plot(1:length(THRESHOLDS),m1_ow1,'bo-','linewidth',1.2);
h11 = plot(1:length(THRESHOLDS),m2_ow1,'ro-','linewidth',1.2);
h12 = plot(1:length(THRESHOLDS),m1_ow2,'bo:','linewidth',1.2,'markersize',MS);
h13 = plot(1:length(THRESHOLDS),m2_ow2,'ro:','linewidth',1.2,'markersize',MS);

plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m1_ow1-se1_ow1;m1_ow1+se1_ow1],'b-','linewidth',1.2);
plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m2_ow1-se2_ow1;m2_ow1+se2_ow1],'r-','linewidth',1.2);
plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m1_ow2-se1_ow2;m1_ow2+se1_ow2],'b:','linewidth',1.2,'markersize',MS);
plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m2_ow2-se2_ow2;m2_ow2+se2_ow2],'r:','linewidth',1.2,'markersize',MS);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
legend([h10 h11 h12 h13],{'Scenes: basic','Reading: basic','Scenes: overweighted','Reading: overweighted'},'box','off','Location','East');
xlim([0.5 12.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);

nbadcomps_exp1_hc1_ow1_filt1  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,1,1,:));
nbadcomps_exp1_hc1_ow1_filt2  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,2,1,:));
nbadcomps_exp1_hc1_ow1_filt3  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,3,1,:));
nbadcomps_exp1_hc1_ow1_filt4  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,4,1,:));
nbadcomps_exp1_hc1_ow1_filt5  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,5,1,:));
nbadcomps_exp1_hc1_ow1_filt6  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,6,1,:));
nbadcomps_exp1_hc1_ow1_filt7  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,7,1,:));
nbadcomps_exp1_hc1_ow1_filt8  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,8,1,:));
nbadcomps_exp1_hc1_ow1_filt9  = squeeze(nbadcomps(1,SUBJECTS,HC_LEVEL,9,1,:));

subplot(4,2,6); hold on;
title('Number of rejected ICs by filter setting')
plot([11.5 11.5],[0 46],'k:') % dotted border: auto-thresh
MS = 4;
hf1 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt1),'ro-','linewidth',0.8);
hf2 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt2),'bo-','linewidth',0.8);
hf3 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt3),'go-','linewidth',0.8);
hf4 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt4),'mo-','linewidth',0.8);
hf5 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt5),'ro-','linewidth',1.2);
hf6 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt6),'bo-','linewidth',1.2);
hf7 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt7),'go-','linewidth',1.2);
hf8 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt8),'mo-','linewidth',1.2);
hf9 = plot(1:length(THRESHOLDS),mean(nbadcomps_exp1_hc1_ow1_filt9),'co-','linewidth',0.8);

% plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m1_ow1-se1_ow1;m1_ow1+se1_ow1],'b-','linewidth',1.2);
% plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m2_ow1-se2_ow1;m2_ow1+se2_ow1],'r-','linewidth',1.2);
% plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m1_ow2-se1_ow2;m1_ow2+se1_ow2],'b:','linewidth',1.2,'markersize',MS);
% plot([1:length(THRESHOLDS);1:length(THRESHOLDS)],[m2_ow2-se2_ow2;m2_ow2+se2_ow2],'r:','linewidth',1.2,'markersize',MS);

xlabel('Variance ratio threshold used')
ylabel('N rejected ICs');
legend([h10 h11 h12 h13],{'Scenes: basic','Reading: basic','Scenes: overweighted','Reading: overweighted'},'box','off','Location','East');
xlim([0.5 12.5])
ylim([0 NCHANS_EEG-1])
plot([0 12.5],[45 45],'k:')  % max. number of possible ICs
% x-labels
set(gca,'XTick',t2plot_numbs);
set(gca,'XTickLabel',t2plot_names);

clear s m2 m2 h1 h2 hs1 hs2 r1 r2 se1_* s2_* x1 x2

% title for figure
[axt,ht]=suplabel('Variance threshold and rejected ICA components','t');
set(ht,'FontSize',18)

if EXPORTFIG
    % export figure
    set(f3,'Position',[0 0 1200 1000])
    export_fig(f3,'figs/new_corrected_plusmsec/figure_thresholds','-eps','-transparent','-painters')
end

% #################################################################################################################

%% -- SUMMARY: Correction results for different *thresholds* (at filter 2 Hz)
LC_TO_PLOT       = 8            % = 2.0 Hz
LW               = 1.1;
LW2              = 0.8;         % linewidth error bars
THRESHOLDLABELS  = 0.5:0.1:2.1;
XLINE            = 6 % = 1.0    % where to put x-line
THRESH2PLOT   = [1:12];         % 0.5:1.5 and Auto
THRESH2PLOT_X = [1:12];

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
m1_cr_0  = mean(abs(metricsac_exp1_ow1_2Hz),1);
m1_cr_1  = mean(abs(metricsac_exp1_ow2_2Hz),1);
sd1_cr_0 =  std(abs(metricsac_exp1_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd1_cr_1 =  std(abs(metricsac_exp1_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));
m1_sp_0  = mean(abs(metricsp_exp1_ow1_2Hz),1)
m1_sp_1  = mean(abs(metricsp_exp1_ow2_2Hz),1)
sd1_sp_0 =  std(abs(metricsp_exp1_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd1_sp_1 =  std(abs(metricsp_exp1_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));
m1_di_0  = mean(abs(metricstim_exp1_ow1_2Hz),1);
m1_di_1  = mean(abs(metricstim_exp1_ow2_2Hz),1);
sd1_di_0 =  std(abs(metricstim_exp1_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd1_di_1 =  std(abs(metricstim_exp1_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));

% reading
m2_cr_0  = mean(abs(metricsac_exp2_ow1_2Hz),1);
m2_cr_1  = mean(abs(metricsac_exp2_ow2_2Hz),1);
sd2_cr_0 =  std(abs(metricsac_exp2_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd2_cr_1 =  std(abs(metricsac_exp2_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));
m2_sp_0  = mean(abs(metricsp_exp2_ow1_2Hz),1)
m2_sp_1  = mean(abs(metricsp_exp2_ow2_2Hz),1)
sd2_sp_0 =  std(abs(metricsp_exp2_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd2_sp_1 =  std(abs(metricsp_exp2_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));
m2_di_0  = mean(abs(metricstim_exp2_ow1_2Hz),1);
m2_di_1  = mean(abs(metricstim_exp2_ow2_2Hz),1);
sd2_di_0 =  std(abs(metricstim_exp2_ow1_2Hz),0,1) / sqrt(length(SUBJECTS));
sd2_di_1 =  std(abs(metricstim_exp2_ow2_2Hz),0,1) / sqrt(length(SUBJECTS));

% f12 = figure('name','Effect of rejection threshold (at 2 Hz)')

% % scenes
% subplot(2,3,1); hold on; title('Corneoretinal artifact');
% h1 = plot([THRESH2PLOT_X],m1_cr_0,'ko-','linewidth',LW)
% h2 = plot([THRESH2PLOT_X],m1_cr_1,'ro-','linewidth',LW)
% % h1 = plot([THRESH2PLOT_X], mean(metricsac_exp1_ow1_2Hz,1),'bx:','linewidth',LW)
% % h2 = plot([THRESH2PLOT_X], mean(metricsac_exp1_ow2_2Hz,1),'ro-','linewidth',LW)
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_cr_0-sd1_cr_0;m1_cr_0+sd1_cr_0],'k-','linewidth',LW2);
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_cr_1-sd1_cr_1;m1_cr_1+sd1_cr_1],'r-','linewidth',LW2);
% legend([h1 h2],{'basic','overweighted'},'box','off')
% set(gca,'Xtick',THRESH2PLOT_X);
% set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
% xlim([0.5 12.5])
% ylim([0 4])
% xlabel('Threshold')
% ylabel('Frontal lateralization [µV]')
% % axis square
% ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
% 
% subplot(2,3,2); hold on; title('Spike potential');
% h1 = plot([THRESH2PLOT_X],m1_sp_0,'ko-','linewidth',LW)
% h2 = plot([THRESH2PLOT_X],m1_sp_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_sp_0-sd1_sp_0;m1_sp_0+sd1_sp_0],'k-','linewidth',LW2);
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_sp_1-sd1_sp_1;m1_sp_1+sd1_sp_1],'r-','linewidth',LW2);
% set(gca,'Xtick',THRESH2PLOT_X);
% set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
% xlim([0.5 12.5])
% ylim([0 2.1])
% xlabel('Threshold')
% ylabel('GFP [µV]')
% % axis square
% ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(4,2,7); hold on; title('ERP distortion');
h1 = plot([THRESH2PLOT_X],m1_di_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT_X],m1_di_1,'ro-','linewidth',LW)
plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_di_0-sd1_di_0;m1_di_0+sd1_di_0],'k-','linewidth',LW2);
plot([THRESH2PLOT_X;THRESH2PLOT_X],[m1_di_1-sd1_di_1;m1_di_1+sd1_di_1],'r-','linewidth',LW2);set(gca,'Xtick',THRESH2PLOT_X);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([0.5 12.5])
ylim([0 2.1])
xlabel('Variance ratio threshold used')
ylabel('GFP of distortion [µV]')
% axis square
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

% % reading
% subplot(2,3,4); hold on; title('Corneoretinal artifact');
% h1 = plot([THRESH2PLOT_X],m2_cr_0,'ko-','linewidth',LW)
% h2 = plot([THRESH2PLOT_X],m2_cr_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_cr_0-sd2_cr_0;m2_cr_0+sd2_cr_0],'k-','linewidth',LW2);
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_cr_1-sd2_cr_1;m2_cr_1+sd2_cr_1],'r-','linewidth',LW2);
% set(gca,'Xtick',THRESH2PLOT_X);
% set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
% xlim([0.5 12.5])
% ylim([0 4])
% xlabel('Threshold')
% ylabel('Frontal lateralization [µV]')
% % axis square
% ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')
% 
% subplot(2,3,5); hold on; title('Spike potential');
% h1 = plot([THRESH2PLOT_X],m2_sp_0,'ko-','linewidth',LW)
% h2 = plot([THRESH2PLOT_X],m2_sp_1,'ro-','linewidth',LW)
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_sp_0-sd2_sp_0;m2_sp_0+sd2_sp_0],'k-','linewidth',LW2);
% plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_sp_1-sd2_sp_1;m2_sp_1+sd2_sp_1],'r-','linewidth',LW2);
% set(gca,'Xtick',THRESH2PLOT_X);
% set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
% xlim([0.5 12.5])
% ylim([0 2.1])
% xlabel('Threshold')
% ylabel('GFP [µV]')
% % axis square
% ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

subplot(4,2,8); hold on; title('ERP distortion');
h1 = plot([THRESH2PLOT_X],m2_di_0,'ko-','linewidth',LW)
h2 = plot([THRESH2PLOT_X],m2_di_1,'ro-','linewidth',LW)
plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_di_0-sd2_di_0;m2_di_0+sd2_di_0],'k-','linewidth',LW2);
plot([THRESH2PLOT_X;THRESH2PLOT_X],[m2_di_1-sd2_di_1;m2_di_1+sd2_di_1],'r-','linewidth',LW2);set(gca,'Xtick',THRESH2PLOT_X);
set(gca,'XTickLabel',{'0.5','','0.7','','0.9','','1.1','','1.3','','1.5','','1.7','','1.9','','auto'});
xlim([0.5 12.5])
ylim([0 2.1])
xlabel('Variance ratio threshold used')
ylabel('GFP of distortion [µV]')
% axis square
ylimits = ylim; plot([XLINE XLINE],[ylimits(1) ylimits(2)],'k:')

% title for figure
[axt,ht]=suplabel('Correction quality (filter = 2 Hz)','t');
set(ht,'FontSize',18)
%[axty,hty]=suplabel('READING     /     SCENES','y');
% set(hty,'FontSize',18)

if EXPORTFIG
    % export figure
    set(f12,'Position',SCREENSIZE)
    export_fig(f12,'figs/new_corrected_plusmsec/figure_summary_thresholds','-eps','-transparent','-painters')
end