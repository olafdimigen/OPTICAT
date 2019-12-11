%% PLOT RESULTS OF "OPTICA" STUDY
% olaf.dimigen@hu-berlin.de, 2018
clear

load('Y:/OPTICA/results/opticaresults.mat','opticaresults')
addpath('M:/Dropbox/_subfunc_master')
addpath M:/Dropbox/eeglab14_1_1b; eeglab; close

SUBJECTS   = 1:12
NCHANS_EEG = 46

%% get grand-average corrected ERPs
n_scene = [];
n_read  = [];

% pre-allocate
SCENE_SACC_CORR = NaN(50,400,2,20,2);
SCENE_STIM_CORR = NaN(50,200,2,20,2);
READ_SACC_CORR  = NaN(50,400,2,20,2);
READ_STIM_CORR  = NaN(50,200,2,20,2);

SCENE_SACC_ORIG = NaN(50,400,12);
SCENE_STIM_ORIG = NaN(50,200,12);
READ_SACC_ORIG  = NaN(50,400,12);
READ_STIM_ORIG  = NaN(50,200,12);

%% get and average original, uncorrected data across subjects
for s = SUBJECTS
    
    fprintf('\nLoading uncorrected data. Subject: %i',s);
    
    % scenes get original data (no ICA)
    load(sprintf('Y:/OPTICA/scenes/correcteddata_new/correcteddata_subj_%i_hc_2_lc_1_ow_1_thresh_7.mat',s),'EEG_sac_RX','EEG_stim_nosacX');
    SCENE_STIM_ORIG(:,:,s) = EEG_stim_nosacX.data(:,1:200);
    SCENE_SACC_ORIG(:,:,s) = EEG_sac_RX.data;
    clear EEG*
    
    % reading: get original data (no ICA)
    load(sprintf('Y:/OPTICA/reading/correcteddata_new/correcteddata_subj_%i_hc_2_lc_1_ow_1_thresh_7.mat',s),'EEG_sac_RX','EEG_stim_nosacX');
    READ_STIM_ORIG(:,:,s) = EEG_stim_nosacX.data(:,1:200);
    READ_SACC_ORIG(:,:,s) = EEG_sac_RX.data;
    clear EEG*
    
    % count number "stim-nosac" epochs
    n_scene = [n_scene; opticaresults(1).subj(s).hc(1).lc(1).ow(1).nepochs_stim_nosac];
    n_read  = [n_read;  opticaresults(2).subj(s).hc(1).lc(1).ow(1).nepochs_stim_nosac];
end
% average over subjects
SCENE_SACC_ORIG  = squeeze(mean(SCENE_SACC_ORIG,3));
SCENE_STIM_ORIG  = squeeze(mean(SCENE_STIM_ORIG,3));
READ_SACC_ORIG   = squeeze(mean(READ_SACC_ORIG ,3));
READ_STIM_ORIG   = squeeze(mean(READ_STIM_ORIG ,3));

fprintf('\nStim-nosac epochs left per subject:\n')
[min(n_scene) max(n_scene) mean(n_scene) std(n_scene)]
[min(n_read)  max(n_read)  mean(n_read)  std(n_read) ]

%% get and average ICA-corrected data
for hc_level = 1:2
    
    for lc_level = 1:20
        
        for ow_level = 1:2
            
            for s = SUBJECTS
                fprintf('\nHC: %i, LC: %i, OW: %i, Subj: %i',hc_level, lc_level, ow_level ,s);
                
                % scene dataset: load subject ERPs
                load(sprintf('Y:/OPTICA/scenes/correcteddata_new/correcteddata_subj_%i_hc_%i_lc_%i_ow_%i_thresh_7.mat',s,hc_level,lc_level,ow_level),'EEG_sac_R2','EEG_stim_nosac2');
                
                tmp_m1_sacc(:,:,s) = EEG_sac_R2.data;
                tmp_m1_stim(:,:,s) = EEG_stim_nosac2.data(:,1:200);
                clear EEG_*
                
                % reading dataset: load subject ERPs
                load(sprintf('Y:/OPTICA/reading/correcteddata_new/correcteddata_subj_%i_hc_%i_lc_%i_ow_%i_thresh_7.mat',s,hc_level,lc_level,ow_level),'EEG_sac_R2','EEG_stim_nosac2');
                
                tmp_m2_sacc(:,:,s) = EEG_sac_R2.data;
                tmp_m2_stim(:,:,s) = EEG_stim_nosac2.data(:,1:200);
                clear EEG_*
                
                % Now average across subjects
                % so only necessary to store grand average for plotting
                % otherwise, a 50x400x12x9x3 matrix has to be maintained in RAM
                if s == 12
                    SCENE_SACC_CORR(:,:,hc_level,lc_level,ow_level) = mean(tmp_m1_sacc,3);
                    SCENE_STIM_CORR(:,:,hc_level,lc_level,ow_level) = mean(tmp_m1_stim,3);
                    
                    READ_SACC_CORR(:,:,hc_level,lc_level,ow_level)  = mean(tmp_m2_sacc,3);
                    READ_STIM_CORR(:,:,hc_level,lc_level,ow_level)  = mean(tmp_m2_stim,3);
                end
                
                if hc_level == 1 & lc_level == 1 & ow_level == 2
                    %% MSEC
                    load(sprintf('Y:/OPTICA/scenes/correcteddata_new/MSEC/correcteddata_subj_%i_MSEC.mat',s),'EEG_sac_RM','EEG_stim_nosacM'); % MSEC-corrected data
                    
                    SCENE_SACC_CORR_MSEC(:,:,s) = EEG_sac_RM.data;
                    SCENE_STIM_CORR_MSEC(:,:,s) = EEG_stim_nosacM.data(:,1:200);
                    clear EEG_sac_RM EEG_stim_nosacM
                    
                    load(sprintf('Y:/OPTICA/reading/correcteddata_new/MSEC/correcteddata_subj_%i_MSEC',s),'EEG_sac_RM','EEG_stim_nosacM'); % MSEC-corrected data
                    
                    READ_SACC_CORR_MSEC(:,:,s) = EEG_sac_RM.data;
                    READ_STIM_CORR_MSEC(:,:,s) = EEG_stim_nosacM.data(:,1:200);
                    clear EEG_sac_RM EEG_stim_nosacM
                    
                end
            end % subj
        end % ow_level
    end % lc level
end % hc level

%% MSEC
% SCENE_SACC_CORR_MSEC = squeeze(mean(SCENE_SACC_CORR_MSEC,3));
% SCENE_STIM_CORR_MSEC = squeeze(mean(SCENE_STIM_CORR_MSEC,3));
% READ_SACC_CORR_MSEC  = squeeze(mean(READ_SACC_CORR_MSEC,3));
% READ_STIM_CORR_MSEC  = squeeze(mean(READ_STIM_CORR_MSEC,3));

whos SCENE* READ*

%% save corrected and uncorrected (orig.) ERPs
save Y:/OPTICA/results/results_for_plotting.mat SCENE* READ*