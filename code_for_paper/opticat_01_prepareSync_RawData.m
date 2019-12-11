%% Prepare the raw data and make it directly comparable between reading and scene viewing exp.
% -- synchronize the EEG data with the eye-tracking data
% -- select a subset of 46 electrodes from dataset 2 (reading) to make it comparable to dataset 1 with 46 electrodes
% -- filter both datasets at ~0.016 Hz (~ 10 s TC) to make them comparable (reading data was aquired as DC)

clear; close all;

addpath M:/Dropbox/_subfunc_master % include special filter (DC corr) that ignores "boundary" events from amplifier resets (DC corrections) 
addpath M:/Dropbox/eeglab14_1_1b
addpath(genpath('M:\Dropbox\erplab7.0.0\'))

syncres = []; % store synch quality results

%% dataset loop: got tru experiements
% 1: scene viewing
% 2: reading
for dataset = 1:2
    
    switch dataset
        case 1 % scenes
            path_raw       = 'Y:/OPTICA/scenes/raw/';
            path_et        = 'Y:/OPTICA/scenes/eyetrack/';
            path_rawset_et = 'Y:/OPTICA/scenes/raw_set_et/';
            SUBJECTS       = [1:12];
            
        case 2 % reading
            path_raw       = 'Y:/OPTICA/reading/raw/';
            path_et        = 'Y:/OPTICA/reading/eyetrack/';
            path_rawset_et = 'Y:/OPTICA/reading/raw_set_et/';
            SUBJECTS       = [1:12];
    end
    
    % subject loop
    for s = SUBJECTS
        
        fprintf('\n\n\nWorking on dataset %i, subject %i...\n',dataset,s);
        
        %% parse eye track
        switch dataset
            
            case 1 % scenes
                eegfilename = sprintf('MSDARK5_%02i.vhdr',s);
                etfilename  = sprintf('%sMSDARK5_%02i Samples.txt',path_et,s);
                etparsed    = sprintf('%sMSDARK5_%02i Samples.mat',path_et,s);
                start_sync_trigger = 103;
                
                % exception for subj 11, different sync triggers
                if ismember(s,[2 11])
                    end_sync_trigger = 203;
                else
                    end_sync_trigger = 99;
                end
                
            case 2 % reading
                eegfilename = sprintf('%02i.vhdr',s);
                etfilename  = sprintf('%sEL5_%i Samples.txt',path_et,s);
                etparsed    = sprintf('%sEL5_%i Samples.mat',path_et,s);
                start_sync_trigger = 200;  end_sync_trigger = 202;
                
        end
        
        %% load original EEG recording
        EEG = pop_loadbv(path_raw,eegfilename);
        
        %% remove channels for reading data (originally 63+1 chans)
        if dataset == 2
            % Reduce reading data to only 46 channels (like the scene data)
            % The montage is overlapping except for 4 electrodes
            KEEP_CHANS = [01 02 03 04 05 06 07 08 10 12 13 15 17 19 21 22 23 24 26 27 28 29 31 33 35 37 38 39 41 42 43 44 46 48 50 52 53 54 56 58 59 60 61 62 63];
            EEG = pop_select(EEG,'channel', KEEP_CHANS);
        end
        
        % low-pass filter all raw data once at 100 Hz
        EEG = pop_eegfiltnew(EEG,[],100);
                
        %% do filtering at 0.0159 Hz (10 sec time constant) to make both datasets comparable
        % (user filter that filters outside of DC boundaries for reading data)
        % the following calls the standard ERPLAB butterworth filter
        % but ignores data breaks (boundary events) resulting from
        % occasional DC corrections of the DC amplifier (reading data was aquired as DC data)
        % that would otherwise potentially cause filtering artifacts
        % Butterworth filter was used here because pop_eegfiltnew too slow to HP-filter at such a low cutoff
        LOWCUTOFF = 0.015915494;
        EEG  = pop_basicfilter_ignore_dccorr(EEG,LOWCUTOFF,400); % call ERPLAB butterworth filter
        
        %% create "zero" channel A1 (for later re-referencing to average ref)
        % but don't include this channel in HP-filtering and in ICAs
        % (it is only added to make average-referencing easier later)
        EEG.nbchan = 46;
        EEG.data(46,:) = 0;
        EEG.chanlocs(46).labels = 'A1';
        EEG = eeg_checkset(EEG);
        
        %% add channel locations
        EEG = pop_chanedit(EEG,'lookup','M:\\Dropbox\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
         
        %% parse eye track (already done)
        if dataset == 1 & s == 2
           et = parsesmi(etfilename,etparsed,'MYKEYWORD'); % for one subject (2), data was synchronized via "keywords" in ET (instead of shared triggers)
        else
           et = parsesmi(etfilename,etparsed); % already done
        end
        
        %% synchronize ET data & add ET channels
        % using EYE-EEG functions
        switch dataset
            case 1               
                EEG = pop_importeyetracker(EEG, etparsed,[start_sync_trigger end_sync_trigger],[7:10],{'Lx' 'Ly' 'Rx' 'Ry'},0,1,0,1);
                close % sync control figure
            case 2
                EEG = pop_importeyetracker(EEG, etparsed,[start_sync_trigger end_sync_trigger],[3:6],{'Lx' 'Ly' 'Rx' 'Ry'},0,1,0,1);
                close % sync control figure
        end
                
        % remember quality of ET/EEG synchonization
        syncres = [syncres, EEG.etc.eyetracker_syncquality(:,2)];
        
        %% NEW, MAY 2018
        % add markers for "bad ET" invervals with out-of-range ET data
        % considered during saccade detection
        ET_CHAN_L = [47:48];  
        ET_CHAN_R = [49:50];    
        TYPE      = 2;        % add "BAD_ET"  markers
        SMP_EXCLUDED_AROUND_MISSING = 25; % 50 ms, ET recording in reading exp. only starts shortly before each stim onset
        EEG = pop_rej_eyecontin(EEG,[ET_CHAN_L ET_CHAN_R] ,[1 1 1 1] ,[1200 800 1200 800],SMP_EXCLUDED_AROUND_MISSING,TYPE); % reject data where gaze position is far beyond screen

        %% detect eye movments 
        % saccade detection parameters (Engbert & Mergenthaler)
        DEG_PER_PIXEL = 0.036;
        VELTHRESH   =   5; % velocity threshold (EXTRA SENSITIVE...)
        MINDUR_SMP  =   5; % minimum saccade duration (in smp)
        PLOTFIG_EM  =   0; % plot EM properties?      
        EEG = pop_detecteyemovements(EEG,[ET_CHAN_L],[ET_CHAN_R],VELTHRESH,MINDUR_SMP,DEG_PER_PIXEL,1,0,25,2,PLOTFIG_EM,1,1);
        EEG = eeg_checkset(EEG);
                
        %% save copy of this original dataset in EEGLAB format
        savefilename  = sprintf('eeg_%02i.set',s);
        pop_saveset(EEG,'filepath',path_rawset_et,'filename',savefilename);
        
        clear EEG start_* end_* etpars*
        
    end % subject loop
end % dataset loop

% write textfile with sync error
%syncres = [syncres, mean(syncres) std(syncres)];   % add mean and std sync quality
save syncquality.mat syncres

fprintf('\n\nDone.')