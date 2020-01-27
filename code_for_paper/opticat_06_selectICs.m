%% Select ICs, apply correction, evaluate results of different ICA variants
% olaf.dimigen@hu-berlin.de, 2018, update 2019

clear
addpath M:/Dropbox/_subfunc_master/deleteoutliers_FileExchange % outlier detection Grubbs
addpath M:/Dropbox/eeglab14_1_2b; eeglab; close

a = tic; % timer

LOWCUTOFFS  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30]; % passband edges
HICUTOFFS   = [40 100]; % passband edges
THRESHOLDS  = [0.5:0.1:1.6]; % thresholds, 1.6 = outlier-detection

% Note: the "THRESHOLD" of 1.6 means that a method for outlier detection is
% used to identify ocular ICs. This could work well with high-density montages.
% However, I did not pursue this road further, because
% the distribution of variance ratios is not Gaussian (as assumed by the
% method used here). It would be cool to look at this possibility
% (individual thresholds based on outlier detection) more thoroughly in future
% research. TL;DR: In the scripts, the 1.6 is the relict of an early exploration.
MS_WITHOUT_SAC_SMP = 100;  % here: in samples (not ms)
OUTLIER_PVALUE     = 0.10; % p-value for Grubbs outlier detection method (not pursued further, assumes Gaussian distrib.)

SUBJECTS  = 1:12

%% load ICA results
load('Z:/OPTICA/results/icaresults.mat', 'icaresults')

% dataset loop
for dataset = 1:2
    
    switch dataset
        case 1 % scenes
            path = 'Z:/OPTICA/scenes/benchmarkdata_new/final_epochs';
        case 2 % reading
            path = 'Z:/OPTICA/reading/benchmarkdata_new/final_epochs';
    end
    
    % subject loop
    for s = SUBJECTS
        
        subtime = tic;
        
        fprintf('\n\n\n');
        warning('Subject %i!',s)
        
        %% load benchmark data: EEG_stim, EEG_stim_nosac, EEG_sac_L, EEG_sac_R, EEG_sac_U, EEG_sac_D
        filename = sprintf('%s/benchmarkdata_subj_%i.mat',path,s);
        load(filename,'EEG_stim','EEG_stim_nosac','EEG_sac_R');
        NCHANS_EEG = 45;
        
        %% now:
        % 1. apply ICA weights to the original, uncorrected, unfiltered "benchmark" data
        % 2. apply different thresholds for IC rejection
        % 3. extract objective quality criteria for stim & sac-locked epochs
        
        %% load current ICA weights (LP filter x HP filter x "ow" (basic/extended of data)
        hc_level = 1;
        for hc_freq = HICUTOFFS
            
            lc_level = 1;
            for lc_freq = LOWCUTOFFS
                
                for ow_level = 1:2
                    
                    fprintf('\n\nGetting weights for dataset: %i, subj: %i, High cutoff: %.1f, Low cutoff: %.1f, overweighting: %i...',dataset,s,HICUTOFFS(hc_level),LOWCUTOFFS(lc_level),ow_level);
                    
                    %% get weights matrix + sphering matrix
                    wts = allicaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).wts;
                    sph = allicaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).sph;
                    
                    
                    %% IMPORTANT! Imported ICA info needs to be deleted if present
                    EEG_stim.icaact            = [];
                    EEG_stim.icasphere         = [];
                    EEG_stim.icaweights        = [];
                    EEG_stim.icachansind       = [];
                    EEG_stim.icawinv           = [];
                    EEG_sac_R.icaact           = [];
                    EEG_sac_R.icasphere        = [];
                    EEG_sac_R.icaweights       = [];
                    EEG_sac_R.icachansind      = [];
                    EEG_sac_R.icawinv          = [];
                    EEG_stim_nosac.icaact      = [];
                    EEG_stim_nosac.icasphere   = [];
                    EEG_stim_nosac.icaweights  = [];
                    EEG_stim_nosac.icachansind = [];
                    EEG_stim_nosac.icawinv     = [];
                    
                    % put in new ICA weights & sphering matrix
                    EEG_stim.icasphere   = sph;
                    EEG_stim.icaweights  = wts;
                    EEG_stim.icachansind = 1:NCHANS_EEG;
                    
                    EEG_stim_nosac.icasphere   = sph;
                    EEG_stim_nosac.icaweights  = wts;
                    EEG_stim_nosac.icachansind = 1:NCHANS_EEG;
                    
                    EEG_sac_R.icasphere   = sph;
                    EEG_sac_R.icaweights  = wts;
                    EEG_sac_R.icachansind = 1:NCHANS_EEG;
                    
                    % let EEGLAB re-compute EEG.icaact & EEG.icawinv
                    EEG_stim       = eeg_checkset(EEG_stim);
                    EEG_sac_R      = eeg_checkset(EEG_sac_R);
                    EEG_stim_nosac = eeg_checkset(EEG_stim_nosac);
                    clear sph wts
                    
                    %% get variance-ratio-table once (do not change data)
                    plotfig        = false;
                    icplotmode     = 4;
                    dummythreshold = 1.0; % used just to get variance ratio table
                    
                    if exist('varratiotable')
                        clear varratiotable
                    end
                    
                    % get variance ratio table
                    % (use saccade window from -10 ms before sacc until sacc offset)
                    [~, varratiotable50] = pop_eyetrackerica(EEG_stim,'saccade','fixation',[5 0], dummythreshold,3,plotfig,icplotmode);
                    
                    % save variance ratio table
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).varratio50 = varratiotable50(:,3);
                    
                    % continue working with window: -10 to 0 ms
                    varratiotable = varratiotable50;
                    
                    %% get threshold cumulative *function*
                    varratio = varratiotable(:,3);
                    ratios = 0:0.1:10;
                    for r = 1:length(ratios)
                        n_rejected(r) = sum(varratio > ratios(r));
                    end
                    
                    %% save sac/fix variances and varratios for this dataset/subj/lc_freq/version
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).varratiotable      = varratiotable(:,1:2); % column 3 is varratio
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).varratio           = varratio;
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).n_rejected         = n_rejected;
                    % add no. of epochs in TEST dataset
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).nepochs_stim       = size(EEG_stim.data,3);
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).nepochs_stim_nosac = size(EEG_stim_nosac.data,3);
                    opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).nepochs_sacR       = size(EEG_sac_R.data,3);
                    
                    % UNCORRECTED DATA: keep only AVERAGES from here on
                    % make copy called "X"
                    EEG_stimX.data         = mean(EEG_stim.data,3);
                    EEG_stim_nosacX.data   = mean(EEG_stim_nosac.data,3);
                    EEG_sac_RX.data        = mean(EEG_sac_R.data,3);
                    EEG_stimX.trials       = 1;
                    EEG_stim_nosacX.trials = 1;
                    EEG_sac_RX.trials      = 1;
                    
                    %% apply different rejection thresholds to remove different ICs
                    thresholdlevel = 1;
                    for thresh = THRESHOLDS
                        
                        %fprintf('\nApplying var. threshold of: %.1f...', thresh);
                        if thresh <= 1.5
                            
                            %% find ICs whos ratio exceeds current threshold
                            badcomps = varratio > thresh;
                            badcomps = badcomps';
                            
                        elseif thresh == 1.6 % not reported/used
                            
                            % special code: detect ocular ICs as outliers
                            % in terms of variance ratio
                            % not pursued further and not reported in paper
                            % and probably only makes sense with many, many
                            % EEG channels/ICs
                            
                            %% detect bad ICs as "outliers" = adaptive method
                            [~,ix_badcomps,badvalues] = deleteoutliers(varratio,OUTLIER_PVALUE);
                            
                            % delete badcomps with varratio < median (=outliers to wrong side!)
                            ix_smallval = varratio(ix_badcomps) < median(varratio);
                            ix_badcomps(ix_smallval) = [];
                            
                            % remember which threshold this "adaptive" outlier algorithm actually corresponded to
                            badIC_withMinRatio = min(varratio(ix_badcomps));
                            opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).smallestOutlier = badIC_withMinRatio;
                            
                            % convert to bool vector
                            badcomps = zeros(1,length(varratio));
                            badcomps(ix_badcomps) = 1;
                        else
                            error('threshold definition problem')
                        end
                        
                        %fprintf('\n------------------------------')
                        %fprintf('\n%i of %i ICs exceed variance ratio of %.3f\n',sum(badcomps),length(badcomps),thresh);
                        
                        %% create copies with removed bad ICs
                        EEG_stim.reject.gcompreject       = badcomps;
                        EEG_stim_nosac.reject.gcompreject = badcomps;
                        EEG_sac_R.reject.gcompreject      = badcomps;
                        
                        if ~all(badcomps) % not all ICs are bad/ocular
                            
                            % actually remove ICs from data
                            EEG_stim2       = pop_subcomp(EEG_stim,find(badcomps),0);
                            EEG_stim2       = eeg_checkset(EEG_stim2);
                            
                            EEG_stim_nosac2 = pop_subcomp(EEG_stim_nosac,find(badcomps),0);
                            EEG_stim_nosac2 = eeg_checkset(EEG_stim_nosac2);
                            
                            EEG_sac_R2      = pop_subcomp(EEG_sac_R,find(badcomps),0);
                            EEG_sac_R2      = eeg_checkset(EEG_sac_R2);
                            
                        else % all ICs are bad/ocular, no output is left (this can happen at very low thresholds): fill output with zeros
                            
                            % EEG_stim
                            nsmp = size(EEG_stim.data,2);
                            ntrl = size(EEG_stim.data,3);
                            EEG_stim2 = EEG_stim;
                            EEG_stim2.data(1:NCHANS_EEG,1:nsmp,1:ntrl)    = zeros(NCHANS_EEG,nsmp,ntrl); % no ICs left so == 0
                            EEG_stim2.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:) = EEG_stim.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:); % keep ET channels
                            
                            % EEG_stim_nosac
                            nsmp = size(EEG_stim_nosac.data,2);
                            ntrl = size(EEG_stim_nosac.data,3);
                            EEG_stim_nosac2 = EEG_stim_nosac;
                            EEG_stim_nosac2.data(1:NCHANS_EEG,1:nsmp,1:ntrl)    = zeros(NCHANS_EEG,nsmp,ntrl); % no ICs left so == 0
                            EEG_stim_nosac2.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:) = EEG_stim_nosac.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:); % keep ET channels
                            
                            % EEG_sac_R
                            nsmp = size(EEG_sac_R.data,2);
                            ntrl = size(EEG_sac_R.data,3);
                            EEG_sac_R2 = EEG_sac_R;
                            EEG_sac_R2.data(1:NCHANS_EEG,1:nsmp,1:ntrl)    = zeros(NCHANS_EEG,nsmp,ntrl); % no ICs left so == 0
                            EEG_sac_R2.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:) = EEG_sac_R.data(NCHANS_EEG+2:NCHANS_EEG+5,:,:); % keep ET channels
                        end
                        
                        % ICA-CORRECTED DATA: keep only AVERAGES from here on
                        % make a copy called "2"
                        EEG_stim2.data         = mean(EEG_stim2.data,3);
                        EEG_stim_nosac2.data   = mean(EEG_stim_nosac2.data,3);
                        EEG_sac_R2.data        = mean(EEG_sac_R2.data,3);
                        EEG_stim2.trials       = 1;
                        EEG_stim_nosac2.trials = 1;
                        EEG_sac_R2.trials      = 1;
                        
                        
                        %% AVERAGE-REFERENCE CORRECTED DATA
                        % EEG = pop_reref( EEG, [],'exclude',[47:50] ); % EEGLAB command
                        
                        % average-reference data "manually":
                        % 3 x original, uncorrected averages
                        m = mean(EEG_stimX.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_stimX.data(1:NCHANS_EEG+1,:) = EEG_stimX.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        m = mean(EEG_stim_nosacX.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_stim_nosacX.data(1:NCHANS_EEG+1,:) = EEG_stim_nosacX.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        m = mean(EEG_sac_RX.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_sac_RX.data(1:NCHANS_EEG+1,:) = EEG_sac_RX.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        
                        % 3 x ICA-corrected averages
                        m = mean(EEG_stim2.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_stim2.data(1:NCHANS_EEG+1,:) = EEG_stim2.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        m = mean(EEG_stim_nosac2.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_stim_nosac2.data(1:NCHANS_EEG+1,:) = EEG_stim_nosac2.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        m = mean(EEG_sac_R2.data(1:NCHANS_EEG+1,:),1);
                        m = repmat(m,NCHANS_EEG+1,1);
                        EEG_sac_R2.data(1:NCHANS_EEG+1,:) = EEG_sac_R2.data(1:NCHANS_EEG+1,:) - m;
                        clear m
                        
                        %% baseline-correct the averages again
                        % because the baseline could have been (slightly) changed by ICA weight application and by average-referencing
                        
                        % stim-locked
                        nframes = size(EEG_stimX.data,2);
                        basevec_smp = 1:50;  % -100 to 0 ms
                        EEG_stimX.data       = rmbase(EEG_stimX.data,nframes,basevec_smp);
                        EEG_stim_nosacX.data = rmbase(EEG_stim_nosacX.data,nframes,basevec_smp);
                        EEG_stim2.data       = rmbase(EEG_stim2.data,nframes,basevec_smp);
                        EEG_stim_nosac2.data = rmbase(EEG_stim_nosac2.data,nframes,basevec_smp);
                        clear nframes basevec_smp
                        
                        % sac-locked
                        nframes = size(EEG_sac_R2.data,2);
                        %basevec_smp = 51:90; % -100 to -20 ms (without SP)
                        basevec_smp  = 76:95; % -50 to -10 ms (without SP)
                        EEG_sac_RX.data      = rmbase(EEG_sac_RX.data,nframes,basevec_smp);
                        EEG_sac_R2.data      = rmbase(EEG_sac_R2.data,nframes,basevec_smp);
                        
                        % get left-vs-right hemisphere frontal electrodes
                        switch dataset
                            case 1 % scenes
                                elec_L = [01 03 05 08 11 12 16 17]; % LO1/2, IO1/2, Fp1/2, AF7/8, F7/8, F3/4, FT9/10 FC5/6
                                elec_R = [02 04 07 10 15 14 21 20];
                            case 2 % reading
                                elec_L = [01 03 05 08 11 12 16 17]; % LO1/2, IO1/2, Fp1/2, AF7/8, F7/8, F3/4, FT9/10 FC5/6
                                elec_R = [02 04 07 10 15 14 21 20];
                        end
                        
                        %% EXTRACT QUALITY METRICS
                        % note: filter ringing from any LP-filtering (see above) can potentially play a role here (for SP)
                        PRESAC_SAMPLES  = abs(EEG_sac_R.xmin)      * EEG_sac_R.srate;      % 100 smp
                        PRESTIM_SAMPLES = abs(EEG_stim_nosac.xmin) * EEG_stim_nosac.srate; %  50 smp
                        
                        % 1. undercorrection: post-saccadic lateralization
                        % take average lateralization from +10 ms to +200 ms
                        d_sac = EEG_sac_R2.data(elec_R,:)-EEG_sac_R2.data(elec_L,:);
                        d_sac_amp = mean(mean(d_sac(:,PRESAC_SAMPLES+6:PRESAC_SAMPLES+100),2),1); % CR-Artifakt = +10 to +200 ms
                        
                        % 2a. undercorrection: mean GFP in fixed SP window
                        % (mean standard deviation across electrodes, across interval)
                        d_sp          = EEG_sac_R2.data(1:NCHANS_EEG,PRESAC_SAMPLES-3:PRESAC_SAMPLES+3); % -6 to +4 ms
                        d_sp_gfp      = std(d_sp(:,:),[],1);
                        d_sp_avggfp = mean(d_sp_gfp);
                        
                        % control: to establish absolute zero estimate for SP,
                        % check how large the GFP is in another interval of
                        % the pre-sacc ERP that is equally long (7 smp) = 14 ms
                        % and equally distant from the baseline-subtraction interval (which was: -50 to -10 ms)
                        d_cntrl          = EEG_sac_R2.data(1:NCHANS_EEG,68:74); % equivalent interval
                        d_cntrl_gfp      = std(d_cntrl(:,:),[],1);
                        d_cntrl_avggfp   = mean(d_cntrl_gfp); % average GFP across the 7 samples
                        
                        % figure; hold on; title('debugging plot')
                        % subplot(2,1,1)
                        % plot(EEG_sac_RX.data(1:46,:)');
                        % plot(EEG_sac_RX.data([3 4 33],:)','r','linewidth',2)
                        % plot(reog,'b','linewidth',2)
                        % subplot(2,1,2)
                        % plot(EEG_sac_R2.data(1:46,:)'); hold on
                        % plot(EEG_sac_R2.data([3 4 33],:)','r','linewidth',2)
                        
                        % 3. overcorrection: stim-ERP distortion: average GFP of difference uncorrected vs. corrected
                        d_stim = EEG_stim_nosacX.data(1:NCHANS_EEG,PRESTIM_SAMPLES:PRESTIM_SAMPLES+MS_WITHOUT_SAC_SMP) - EEG_stim_nosac2.data(1:NCHANS_EEG,PRESTIM_SAMPLES:PRESTIM_SAMPLES+MS_WITHOUT_SAC_SMP);
                        d_stim_gfp = std(d_stim(:,:),[],1);
                        d_stim_avggfp = mean(d_stim_gfp);
                        
                        %% add more results
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_sac       = d_sac_amp;
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_sp        = d_sp_avggfp;
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_sp3       = d_sp_avggfp_roi;
                        
                        % window-based: control window on "other" side of baseline interval (smp 68:74) to see how much GFP is expected
                        % to obtain a zero-point estimate for the SP amplitude (since GFP will rarely fall to zero)
                        % (this is only used as a reference value for plotting in Figure 3)
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_cntrl     = d_cntrl_avggfp;
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_cntrl_roi = d_cntrl_avggfp_roi;
                        
                        %
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).metric_stim      = d_stim_avggfp;
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).topo_metric_sac  = mean(d_sac,2);
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).topo_metric_sp   = mean(d_sp,2);
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).topo_metric_stim = mean(d_stim,2);
                        
                        % add infos: which/how many bad comps for this threshold?
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).badcomps         = find(badcomps);
                        opticaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(thresholdlevel).nbadcomps        = sum(badcomps);
                        
                        %% save full corrected waveforms --> ONLY for some conditions (threshold = 1.1)
                        if THRESHOLDS(thresholdlevel) == 1.1
                            
                            %% save corrected data: EEG_stim_nosac2, EEG_sac_R2,
                            switch dataset
                                case 1
                                    savepath = sprintf('Z:/OPTICA/scenes/correcteddata_new/correcteddata_subj_%i_hc_%i_lc_%i_ow_%i_thresh_%i.mat',s,hc_level,lc_level,ow_level,thresholdlevel);
                                case 2
                                    savepath = sprintf('Z:/OPTICA/reading/correcteddata_new/correcteddata_subj_%i_hc_%i_lc_%i_ow_%i_thresh_%i.mat',s,hc_level,lc_level,ow_level,thresholdlevel);
                            end
                            % save
                            save(savepath,'EEG_sac_RX','EEG_sac_R2','EEG_stimX','EEG_stim2','EEG_stim_nosacX','EEG_stim_nosac2'); %EEG_stim, EEG_stim_nosac, EEG_sacR
                            
                        end
                        
                        clear d_* EEG_stim2 EEG_nosac2 EEG_sac_R2 % varratio badcomps
                        
                        thresholdlevel = thresholdlevel+1;
                    end % thresh loop
                    
                    clear wts sph
                end % ow_level
                
                lc_level = lc_level+1;
            end % high-pass filter
            
            hc_level = hc_level+1;
            
        end % low-pass filter
        
        %% do surrogate MSEC (Berg & Scherg, 1998) correction and extract the same metrics
        % -- load msec correction matrix (generated in "BESA" from text file)
        % -- convert ERPs to average-reference (since correction matrix is for avg-ref data)
        % -- multiply matrix with averaged ERPs
        % -- extract quality metrics
        
        % delete metrics (from ICA)
        clear d_sac_amp d_sp_avggfp sp_gfp d_stim_avggfp d_sac d_sp d_stim
        
        %% 1. read the BESA matrix
        clear besamatrix
        delimiter = '\t';
        startRow = 2;
        % get filename/format of besamatrix
        
        switch dataset
            case 1
                msecfilename = sprintf('Z:/OPTICA/scenes/msecmatrix/c%02d.matrix',s);
            case 2
                msecfilename = sprintf('Z:/OPTICA/reading/msecmatrix/redone2018_c%02d.matrix',s);
        end
        
        formatSpec = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]'; % 46 channels
        fileID = fopen(msecfilename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID); % close text file
        besamatrix = [dataArray{1:end-1}];
        clear msecfilename delimiter startRow formatSpec fileID dataArray
        fprintf('\nSize of BESA MSEC matrix BEFORE: %i x %i\n',size(besamatrix,1),size(besamatrix,2));
        besamatrix = besamatrix(1:46,1:46); % cut away the 3 MSEC "ART" channels, MAJOR BUGFIX JUNE 2nd, 2018
        fprintf('\nSize of BESA MSEC matrix AFTER: %i x %i\n',size(besamatrix,1),size(besamatrix,2));
        % figure; imagesc(besamatrix), colorbar
        
        
        %% get copies of the uncorrected ERPs (stim, stim_nosac, sac_R)
        EEG_stimM.data       = mean(EEG_stim.data,3);
        EEG_stim_nosacM.data = mean(EEG_stim_nosac.data,3);
        EEG_sac_RM.data      = mean(EEG_sac_R.data,3);
        
        
        %% average-reference them (because BESA-matrix computed on avg-ref'ed data)
        m = mean(EEG_stimM.data(1:46,:),1);
        m = repmat(m,46,1);
        EEG_stimM.data(1:46,:) = EEG_stimM.data(1:46,:) - m;
        clear m
        m = mean(EEG_stim_nosacM.data(1:46,:),1);
        m = repmat(m,46,1);
        EEG_stim_nosacM.data(1:46,:) = EEG_stim_nosacM.data(1:46,:) - m;
        clear m
        m = mean(EEG_sac_RM.data(1:46,:),1);
        m = repmat(m,46,1);
        EEG_sac_RM.data(1:46,:) = EEG_sac_RM.data(1:46,:) - m;
        clear m
        
        %% apply BESA matrix by matrix multiplication
        EEG_stimM.data(1:46,:)       = besamatrix * EEG_stimM.data(1:46,:);
        EEG_stim_nosacM.data(1:46,:) = besamatrix * EEG_stim_nosacM.data(1:46,:);
        EEG_sac_RM.data(1:46,:)      = besamatrix * EEG_sac_RM.data(1:46,:);
        
        %% baseline-correct BESA-corrected data again
        % (baseline could have been changed by MSEC and average-referencing, see ICA)
        % stim-locked
        nframes = size(EEG_stimM.data,2);
        basevec_smp = 1:50;  % -100 to 0 ms
        EEG_stimM.data       = rmbase(EEG_stimM.data,nframes,basevec_smp);
        EEG_stim_nosacM.data = rmbase(EEG_stim_nosacM.data,nframes,basevec_smp);
        clear nframes basevec_smp
        % sac-locked
        nframes = size(EEG_sac_RM.data,2);
        basevec_smp = 76:95; % -50 to -10 ms (without SP)
        EEG_sac_RM.data      = rmbase(EEG_sac_RM.data,nframes,basevec_smp);
        
        %% get left-vs-right hemisphere electrodes
        switch dataset
            case 1 % scenes
                elec_L = [01 03 05 08 11 12 16 17]; % LO1/2, IO1/2, Fp1/2, AF7/8, F7/8, F3/4, FT9/10 FC5/6
                elec_R = [02 04 07 10 15 14 21 20];
            case 2 % reading
                elec_L = [01 03 05 08 11 12 16 17]; % LO1/2, IO1/2, Fp1/2, AF7/8, F7/8, F3/4, FT9/10 FC5/6
                elec_R = [02 04 07 10 15 14 21 20];
        end
        
        %% EXTRACT QUALITY METRICS
        % 1. undercorrection: post-saccadic lateralization
        % take average lateralization from +10 ms to +200 ms
        d_sac = EEG_sac_RM.data(elec_R,:)-EEG_sac_RM.data(elec_L,:);
        d_sac_amp = mean(mean(d_sac(:,PRESAC_SAMPLES+6:PRESAC_SAMPLES+100),2),1); % CR-Artifakt = +10 to +200 ms
        
        % 2a. undercorrection: mean GFP in fixed SP window
        % (mean standard deviation across electrodes, across interval)
        d_sp        = EEG_sac_RM.data(1:46,PRESAC_SAMPLES-3:PRESAC_SAMPLES+4); % -6 to +6 ms
        d_sp_gfp    = std(d_sp(:,:),[],1);
        d_sp_avggfp = mean(d_sp_gfp);
        
        % MSEC: control: to establish absolute zero estimate for SP,
        % check how large the GFP is in another interval of
        % the pre-sacc ERP that is equally long (7 smp) = 14 ms
        % and equally distant from baseline-subtraction interval (-50 to -10 ms)
        d_cntrl          = EEG_sac_RM.data(1:NCHANS_EEG,68:74); % equivalent interval
        d_cntrl_gfp      = std(d_cntrl(:,:),[],1);
        d_cntrl_gfp_roi  = std(d_cntrl([elec_L,elec_R],:),[],1); % compute "control" GFP only in frontal ROI
        % average GFP across the 7 samples
        d_cntrl_avggfp     = mean(d_cntrl_gfp);
        d_cntrl_avggfp_roi = mean(d_cntrl_gfp_roi);
        
        % 3. overcorrection: stim-ERP distortion: average GFP of difference uncorrected vs. corrected
        d_stim = EEG_stim_nosacX.data(1:46,PRESTIM_SAMPLES:PRESTIM_SAMPLES+MS_WITHOUT_SAC_SMP) - EEG_stim_nosacM.data(1:46,PRESTIM_SAMPLES:PRESTIM_SAMPLES+MS_WITHOUT_SAC_SMP);
        d_stim_gfp = std(d_stim(:,:),[],1);
        d_stim_avggfp = mean(d_stim_gfp);
        
        %% add MSEC results
        opticaresults(dataset).subj(s).msec.metric_sac       = d_sac_amp;
        opticaresults(dataset).subj(s).msec.metric_sp        = d_sp_avggfp;
        opticaresults(dataset).subj(s).msec.metric_stim      = d_stim_avggfp;
        
        % update: also store GFP "baseline" (absolute zero)
        opticaresults(dataset).subj(s).msec.metric_cntrl     = d_cntrl_avggfp;
        opticaresults(dataset).subj(s).msec.metric_cntrl_roi = d_cntrl_avggfp_roi;
        
        opticaresults(dataset).subj(s).msec.topo_metric_sac  = mean(d_sac,2);
        opticaresults(dataset).subj(s).msec.topo_metric_sp   = mean(d_sp,2);
        opticaresults(dataset).subj(s).msec.topo_metric_stim = mean(d_stim,2);
        
        %% save corrected data: EEG_stim_nosacM, EEG_sac_RM,
        switch dataset
            case 1
                savepathMSEC = sprintf('Z:/OPTICA/scenes/correcteddata_new/MSEC/correcteddata_subj_%i_MSEC.mat',s);
            case 2
                savepathMSEC = sprintf('Z:/OPTICA/reading/correcteddata_new/MSEC/correcteddata_subj_%i_MSEC.mat',s);
        end
        
        % save
        save(savepathMSEC,'EEG_stimM','EEG_sac_RM','EEG_stim_nosacM','besamatrix');
        
        fprintf('\n\n\n\n')
        tt = toc(subtime);
        
        warning('Subject time in minutes:'), disp(tt/60)
        opticaresults(dataset).subj(s)
        
        
    end % sub loop
end % dataset loop

%% save everything
save Z:/OPTICA/results/opticaresults.mat opticaresults LOWCUTOFFS HICUTOFFS THRESHOLDS

total_script_time = toc(a);
fprintf('\nDone.')