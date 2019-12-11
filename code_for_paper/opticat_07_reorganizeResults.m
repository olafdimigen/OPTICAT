%% Re-format results for easier analysis/plotting
% olaf.dimigen@hu-berlin.de, 2017, updated 2018

clear

% load ICA output
load Z:/OPTICA/results/opticaresults.mat opticaresults LOWCUTOFFS HICUTOFFS THRESHOLDS
fprintf('\nData loaded.')

% opticaresults(1:2).subj(1:12).hc(1:2).lc(1:20).ow(1:2)

% varratiotable:      [46x3 double]
% varratio:           [46x1 double]
% n_rejected:         [1x101 double]
% nepochs_stim:       24
% nepochs_stim_nosac: 12
% nepochs_sacR:       43
%
% thresh:             [1x46 struct]

% opticaresults(1:2).subj(1:12).hc(1:2).lc(1:20).ow(1:2).thresh(1:46)
% metric_sac
% metric_sp
% metric_stim
% badcomps
% nbadcomps
% topo_metric_sac
% topo_metric_sp
% topo_metric_stim

% go tru data, reorganize
for exp = 1:2 % experiment: 1 = scenes, 2 = reading
    
    for s = 1:12 % subjects
        
        for hc_level = 1:2 % passband edge low-pass filter: 40 vs 100 Hz
            
            for lc_level = 1:20  % passband edge high-pass filter (20x)
                
                for ow_level = 1:2 % overweighting (no/yes)
                    
                    fprintf('\nExp: %i. Subj: %i. Highcutofflevel: %i, Lowcutofflevel: %i. OW_level: %i\n',exp,s,hc_level,lc_level,ow_level);
                    
                    % infos about sac/fix ratio of ICs
                    exp_varratio(exp,s,hc_level,lc_level,ow_level,:)    = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).varratio;    % variance ratios for all ICs (46 values)
                    exp_n_rejected1(exp,s,hc_level,lc_level,ow_level,:) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).n_rejected;  % rejection function (101 values)
                    
                    % now collect metrics of undercorr. and overcorr.
                    for t = 1:length(THRESHOLDS)
                        fprintf('%i..',t);
                        nbadcomps(exp,s,hc_level,lc_level,ow_level,t)   = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).nbadcomps;   % number of ICs
                        metric_sac(exp,s,hc_level,lc_level,ow_level,t)  = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).metric_sac;  % residual artifact: cr
                        metric_sp(exp,s,hc_level,lc_level,ow_level,t)   = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).metric_sp;   % residual artifact: sp
                        metric_stim(exp,s,hc_level,lc_level,ow_level,t) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).metric_stim; % distortion
                        
                        % new: baseline GFP as absolute zero point
                        metric_baseGFP(exp,s,hc_level,lc_level,ow_level,t) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).metric_cntrl;     % baseline GFP
                        metric_baseGFP_roi(exp,s,hc_level,lc_level,ow_level,t) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).metric_cntrl_roi; % baseline GFP: ROI
                        
                        switch exp
                            case 1
                                exp1_topos_sac  = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_sac; % half-topo: left vs. right difference
                                exp1_topos_sp   = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_sp;  % topo: stim vs. org-stim
                                exp1_topos_stim = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_stim; % topo: sp
                            case 2
                                exp2_topos_sac  = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_sac;
                                exp2_topos_sp   = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_sp;
                                exp2_topos_stim = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(t).topo_metric_stim;
                        end % switch
                        
                    end % thresholds
                end % ow_level
            end % lc_level
        end % hc_level
        
        %% MSEC
        metric_sac_msec(exp,s)   = opticaresults(exp).subj(s).msec.metric_sac;
        metric_sp_msec(exp,s)    = opticaresults(exp).subj(s).msec.metric_sp;
        metric_stim_msec(exp,s)  = opticaresults(exp).subj(s).msec.metric_stim;
        
        % new: baseline GFP as absolute zero point (Figure 3)
        metric_baseGFP(exp,s)     = opticaresults(exp).subj(s).msec.metric_cntrl;     % baseline GFP
        metric_baseGFP_roi(exp,s) = opticaresults(exp).subj(s).msec.metric_cntrl_roi; % baseline GFP: ROI
        
        
        switch exp
            case 1
                exp1_topos_sac_msec  = opticaresults(exp).subj(s).msec.topo_metric_sac;  % half-topo: left vs. right difference
                exp1_topos_sp_msec   = opticaresults(exp).subj(s).msec.topo_metric_sp;   % topo: stim vs. org-stim
                exp1_topos_stim_msec = opticaresults(exp).subj(s).msec.topo_metric_stim; % topo: sp
            case 2
                exp2_topos_sac_msec  = opticaresults(exp).subj(s).msec.topo_metric_sac;  % half-topo: left vs. right difference
                exp2_topos_sp_msec   = opticaresults(exp).subj(s).msec.topo_metric_sp;   % topo: stim vs. org-stim
                exp2_topos_stim_msec = opticaresults(exp).subj(s).msec.topo_metric_stim; % topo: sp
        end % switch
        
    end % s
end % exp

%% save
save Z:/OPTICA/results/metrics.mat metric* nbadcomps exp1_* THR* LOW* exp2_*
fprintf('\nReorganize Results: Done.')