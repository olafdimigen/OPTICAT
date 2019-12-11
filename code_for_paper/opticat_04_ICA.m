%% Compute lots of ICAs for differently processed training data
% Note: On a single machine, this scripts runs for a few days...
clear; close all

addpath M:/Dropbox/_subfunc_master
addpath M:/Dropbox/eeglab14_1_1b
eeglab, close
addpath M:/Dropbox/binica

% filter constants (Hz)
LOWCUTOFFS  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30]; % passband edges
HICUTOFFS   = [40 100];  % passband edges
OWS         = [0.0 1.0]; % overweighting proportions
PLOT_IC     = false;
SUBJECTS    = [1:12];

%% dataset loop
% 1: scenes
% 2: reading
for dataset = 1:2
    switch dataset
        case 1 % scenes
            path_train       = 'Y:/OPTICA/scenes/trainingdata_new/';
            path_ica         = 'Y:/OPTICA/scenes/ica/';
        case 2 % reading
            path_train       = 'Y:/OPTICA/reading/trainingdata_new/';
            path_ica         = 'Y:/OPTICA/reading/ica/';
    end
    
    % subject loop
    for s = SUBJECTS
        
        tic_subjectloop = tic;
        hc_level = 1;
        
        % filter loop
        for hicutoffs = HICUTOFFS
            
            lc_level = 1;
            
            % filter loop
            for lowcutoffs = LOWCUTOFFS
                
                ow_level = 1;
                
                % ow_level: overweighting factor
                for owfactor = OWS
                    
                    
                    warning('LOADING NEXT DATASET...')
                    fprintf('\nLoading dataset: %i subj: %i lowcut: %.1f Hz highcut: %.1f Hz owfactor %.1f\n\n', dataset,s,lowcutoffs,hicutoffs,owfactor)
                    
                    %% load training data
                    loadfilename = sprintf('traineeg_%02i_locut_%.1f_hicut_%i_ow_%.1f.set',s,lowcutoffs,hicutoffs,owfactor);
                    EEG          = pop_loadset('filename',loadfilename,'filepath',path_train);
                    
                    %% run binary ICA
                    tic_ica = tic;
                    [wts,sph] = binica(EEG.data(1:45,:),'extended',1);
                    
                    time_for_ica = toc(tic_ica);
                    fprintf('\n\nMinutes for ICA: %.1f\n',time_for_ica/60);
                    
                    %% save weights & sphering matrix in big struct
                    totalpoints = EEG.pnts * EEG.trials;
                    pointsperweight = totalpoints./45^2
                    fprintf('\n\nData points per computed ICA weight: %.1f\n',pointsperweight);
                    
                    %% add sph & wts also to dataset (optional)
                    EEG.icasphere   = sph;
                    EEG.icaweights  = wts;
                    EEG.icachansind = 1:45;
                    
                    EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact and EEG.icawinv;
                    [ALLEEG EEG] = eeg_store(ALLEEG,EEG,CURRENTSET); % datasets now includes ICA weights
                    eeglab redraw
                    
                    icaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).wts         = wts;
                    icaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).sph         = sph;
                    icaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).wtsinv      = EEG.icawinv;
                    icaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).icapoints   = totalpoints;
                    icaresults(dataset).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).minutes4ica = (time_for_ica/60); % minutes to run ICA
                    
                    %% show some IC topos
                    if PLOT_IC
                        figtit = sprintf('dataset: %i. subj: %i. lowcutoffs: %.2f Hz. overweighting: %i',dataset,s,lowcutoffs,owfactor);
                        pop_topoplot(EEG,0,[1:45],figtit,[5 9],0,'electrodes','on');
                    end
                    
                    %% save actual datasets with ICA weights (this is optional)
                    savefilename = sprintf('ica_%02i_locut_%.1f_hicut_%i_vers_%i.set',s,lowcutoff,highcutoff,ow);
                    EEG = pop_saveset(EEG,'filename',savefilename,'filepath',path_ica);
                    
                    ow_level = ow_level +1;
                    
                end % ow_level
                
                lc_level = lc_level+1;
            end % lowcutoffs loop
            
            hc_level = hc_level+1;
        end % highcutoffs loop
        
        time_for_subject = toc(tic_subjectloop)
       
    end % subj loop
end % dataset loop

% save again, just to be sure
save Y:/OPTICA/results/icaresults.mat icaresults