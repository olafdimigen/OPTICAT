%% create different versions of training data for ICA
% -- different high-pass filters
% -- different low-pass filters
% -- versions without (1) or with (2+3) appended copies of sac-epochs containing the SP
% olaf.dimigen@hu-berlin.de, 2017

clear 
addpath M:/Dropbox/eeglab14_1_1b
eeglab, close

% filter constants (Hz)
LOWCUTOFFS  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30]; % passband edges
HICUTOFFS   = [40 100];  % passband edges

% length of stimulus-locked epochs
BEFORETRIG  = -0.2;
AFTERTRIG   =  2.8; 

% length of (very short) saccade-locked epochs capturing the SP
BEFORESAC   = -0.020;
AFTERSAC    =  0.010;

% other constants
FREAKTHRESH = 500; % im µV: amplitude threshold for removing/pruning epochs with "freak" artifacts
NCHANS_EEG  =  46; % here including ref channel A1 (contains only zeros at this point)

% ICA data length settings
POINTS_PER_WEIGHT      = 80   % points per weight
OVERWEIGHT_PROPORTION1 = 0.2  % proportion of extra sac-locked data
OVERWEIGHT_PROPORTION2 = 1.0  % proportion of extra sac-locked data
NPOINTS                = 45^2 * POINTS_PER_WEIGHT
OTHER_OW_PROPORTIONS   = [0.1:0.1:0.9] % 1.5 2.0 5.0]
REMOVE_MEAN            = 1;   % remove epoch mean of overweighted SP epochs?

SUBJECTS = 1:12

%% dataset loop
for dataset =  1:2
    
    switch dataset
        case 1 % scenes
            path_filt        = 'Y:/OPTICA/scenes/filt_new/';
            path_train       = 'Y:/OPTICA/scenes/trainingdata_new/';
            trigger          = {'S 11' 'S 12' 'S 13' 'S  1'}; % stim onset/search target onset
            path_allepochs   = 'Y:/OPTICA/scenes/benchmarkdata_new/';
        case 2 % reading
            path_filt        = 'Y:/OPTICA/reading/filt_new/';
            path_train       = 'Y:/OPTICA/reading/trainingdata_new/';
            trigger          = {'S101' 'S102' 'S103' 'S104' 'S151' 'S152' 'S153' 'S154'}; % sentence onset
            path_allepochs   = 'Y:/OPTICA/reading/benchmarkdata_new/';
    end
    
    % subject loop
    for s = SUBJECTS
        
        tic_subloop = tic;
        
        % filter loop
        for hicutoff = HICUTOFFS
            
            % filter loop
            for lowcutoff = LOWCUTOFFS
                
                fprintf('\n\n------------------------------------------------')
                fprintf('\nDataset: %i. Subject %i. Low: %.1f Hz  High: %.1f Hz',dataset,s,lowcutoff,hicutoff);
                fprintf('\n--------------------------------------------------\n\n')
                
                %% load
                loadfilename = sprintf('eeg_%02i_lowcut_%.1fHz_hicut_%iHz.set',s,lowcutoff,hicutoff);
                
                EEG = pop_loadset('filename',loadfilename,'filepath',path_filt);
                [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,1);
                
                chanlocs = EEG.chanlocs; % store these in order to add back later
                
                %% 1. cut epochs around stimuli
                EEG = pop_epoch(EEG,trigger,[BEFORETRIG AFTERTRIG], 'newname', 'Stimlocked_Training_Epochs','epochinfo','yes');
                EEG = applytochannels(EEG,[1:NCHANS_EEG],'pop_rmbase( EEG,[]);'); % EYE-EEG wrapper function, baseline-correct EEG channels across whole epoch
                
                %% remove epochs with "freak" (extreme) voltage outlier values
                % this is the only type of "pruning" we will do
                REJECTSTATUS = 1;
                [EEG,ix] = pop_eegthresh(EEG,1,[1:NCHANS_EEG],-FREAKTHRESH,FREAKTHRESH,BEFORETRIG,AFTERTRIG,1,REJECTSTATUS);
                
                %% save epoched (training) data before length is reduced to k = 80 samples per weight
                savefilename  = sprintf('traineeg_%02i_locut_%.1f_hicut_%i_allepochs.set',s,lowcutoff,hicutoff);
                EEG = pop_saveset(EEG, 'filename',savefilename,'filepath',path_allepochs);  % stim only, all epochs
                
                % -----------------------------------------------------
                % Remove extra data that is longer than the intended length of
                % the training data. This step may seem a bit confusing here,
                % but it's necessary to ensure that extra sacc-locked data that
                % is in the "extended" (=overweighted) version of the training data is really
                % entirely redundant with already existing samples in the non-overweighted data
                maxepoch2keep = ceil(NPOINTS / size(EEG.data,2));
                EEG = pop_select(EEG,'trial',1:maxepoch2keep);
                                
                %% save version of training data without overweighting
                ow_factor = 0.0;
                savefilename  = sprintf('traineeg_%02i_locut_%.1f_hicut_%i_ow_%.1f.set',s,lowcutoff,hicutoff,ow_factor);
                EEG = pop_saveset(EEG, 'filename',savefilename,'filepath',path_train);  % stim only
                                                
                %% create and save training data with overweighting
                ow_factor = 1.0;
                [EEG2] = overweightevents(EEG,'saccade',[BEFORESAC AFTERSAC],ow_factor,REMOVE_MEAN);
                
                %% save training data with overweighting
                savefilename_ow = sprintf('traineeg_%02i_locut_%.1f_hicut_%i_ow_%.1f.set',s,lowcutoff,hicutoff,ow_factor);
                EEG2 = pop_saveset(EEG2,'filename',savefilename_ow, 'filepath',path_train);  % stim + 10% sac
                
                %% explore some other overweighting factors
                % supplementary control analysis
                % at HP-filter settings of 0.1 Hz and 2 Hz try OTHER, intermediate and even more extreme overweighting proportions
                if ismember(lowcutoff,[0.1 2] && hicutoff == 100
                    for owx = OTHER_OW_PROPORTIONS
                        fprintf('\nCreating data with different overweighting proportions...%.1f',owx);
                        [EEG3] = overweightevents(EEG,'saccade',[BEFORESAC AFTERSAC],owx,REMOVE_MEAN);
                        % save
                        savefilename_owx = sprintf('traineeg_%02i_locut_%.1f_hicut_%i_ow_%.1f.set',s,lowcutoff,hicutoff,owx)
                        EEG3 = pop_saveset(EEG3,'filename',savefilename_owx, 'filepath',path_train);  % stim + 10% sac
                    end
                end
            end % highpassfilter
        end % lowpassfilter
        
        time_for_subject = toc(tic_subloop);
        fprintf('\n\n----------------------------------------------------')
        fprintf('\nTime for subject %i: %.1f s',s,time_for_subject);
        fprintf('\n------------------------------------------------------')
        
    end % subject
end % dataset

fprintf('\n\nDone.');