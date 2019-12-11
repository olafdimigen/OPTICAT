%% Prepare benchmark data to test the quality of ocular correction
% olaf.dimigen@hu-berlin.de, 2017
clear

addpath M:/Dropbox/eeglab14_1_1b
eeglab, close

MS_WITHOUT_SAC      =  200;   % required interval without (micro)saccade in stim_nosacc epochs
DEG_PER_PIXEL       =  0.036; % degree of visual angle per pixel, screen distance of 60 cm
GAZETHRESHOLD_PIXEL =  0.5 / DEG_PER_PIXEL;
MAXGAZEJUMP         =  0.2 / DEG_PER_PIXEL;
NCHANS_EEG          =  45;

PLOT_SAC_ANGLES     = 1; % control plot

% dataset 1: scenes
singletrialgazeLX1  = []; % left eye: x
singletrialgazeLY1  = []; % left eye: y
singletrialgazeRX1  = []; % right eye: x
singletrialgazeRY1  = []; % right eye: y
% dataset 2: reading
singletrialgazeLX2  = []; 
singletrialgazeLY2  = [];
singletrialgazeRX2  = [];
singletrialgazeRY2  = [];
counter             =  0;

SUBJECTS = 1:12

% dataset loop
for dataset = 1:2
    
    % load epochated data whose first 80 * 45^2 = 162,000 samples were used for ICA. 
    % Epochs with voltages > 500 µV were removed. Contains saccade events, detected in contin. dataset
    switch dataset
        case 1 % scenes
            path_rawset_et   = 'Y:/OPTICA/scenes/benchmarkdata_new/';
            trigger          = {'S 11' 'S 12' 'S 13'}; % without 'S  1', dont use those
        case 2 % reading
            path_rawset_et   = 'Y:/OPTICA/reading/benchmarkdata_new/';
            trigger          = {'S101' 'S102' 'S103' 'S104' 'S151' 'S152' 'S153' 'S154'}; % includes 8 practice trials = 152 trials
    end
    
    % subject loop
    for s = SUBJECTS
        
        %% load raw EEG data (.set)
        loadfilename_rawset = sprintf('traineeg_%02i_locut_0.0_hicut_100_allepochs.set',s);
        EEG = pop_loadset('filename',loadfilename_rawset,'filepath',path_rawset_et);
        EEG = eeg_checkset(EEG);
        
        %% cut baselined epochs around stimulus onsets
        switch dataset
            case 1
                % scenes
                EEG = pop_epoch( EEG, trigger, [-0.1 2.8], 'newname', 'S-ON epochs', 'epochinfo', 'yes');
                EEG = applytochannels(EEG,[1:NCHANS_EEG],' pop_rmbase( EEG, [-100 0]);'); % 100 ms pre-stim baseline
            case 2
                % reading
                EEG = pop_epoch( EEG, trigger, [-0.1 2.8], 'newname', 'S-ON epochs', 'epochinfo', 'yes');
                EEG = applytochannels(EEG,[1:NCHANS_EEG],' pop_rmbase( EEG, [-100 0]);'); % 100 ms pre-stim baseline
        end
               
        % % reject epochs with extreme EEG (e.g. DC corrections)
        % EEG = pop_eegthresh(EEG,1,[1:NCHANS_EEG],[-500],[500],EEG.times(1)/1000,EEG.times(end)/1000,0,REJECTNOW);
        % EEG = eeg_checkset(EEG);
               
        %% save fixation/saccade properties
        ix = find(ismember({EEG.event.type},{'saccade','fixation'}));
        eyemovements(dataset).subject(s).e = EEG.event(ix);
        
        %% now we need two more epochations:
        % -- around saccades (= SRP), left vs. right hemisphere or sac dir
        % -- around stimulus onsets, without ANY (micro)saccade in first X ms
        % Examplary numbers of observations:
                
        %% 1. Check 1: around stimulus onsets without a sacc within "MS_WITHOUT_SAC"
        % search for first sacc in epoch
        EEG_stim = EEG;
        [sac_lat_1stSaccade, ~] = eeg_getepochevent(EEG_stim,'saccade',[],'latency'); % note: produces output in ms!
        ix_earlysacc = find(sac_lat_1stSaccade < MS_WITHOUT_SAC); % find epochs without saccades
        
        %% EEG_stim_nosac: plot single-trial ET data and onset of first detected saccade
        %         time = -100:2:4998
        %         for t = 1:size(EEG_stim.data,3)
        %             t
        %             figure; hold on;
        %             ylim([200 900]); xlim([-100 400])
        %
        %             plot(time, EEG_stim.data(47:49,:,t));
        %             plot([sac_lat_1stSaccade(t) sac_lat_1stSaccade(t)],[200 900],'r:')
        %             sac_lat_1stSaccade(t)
        %             disp('reject: ')
        %             sac_lat_1stSaccade(t) < MS_WITHOUT_SAC
        %             plot([250 250],[200 900],'k--')
        %             pause
        %             close
        %         end
        
        % *REMOVE* trials with an early saccade event (flag: 'notrial')
        EEG_stim_nosac = pop_select(EEG_stim,'notrial',[ix_earlysacc]);
        EEG_stim_nosac = eeg_checkset(EEG_stim_nosac);
                
        %% Gaze check 2: check for max. gaze deviation within the epoch
        start = 1 % abs(EEG_stim_nosac.times(1)/ (1000/EEG_stim_nosac.srate)) + 1;
        stop  = abs(EEG_stim_nosac.times(1)/ (1000/EEG_stim_nosac.srate)) + (MS_WITHOUT_SAC / (1000/EEG_stim_nosac.srate));
        
        killme = [];
        for t = 1:size(EEG_stim_nosac.data,3)
            
            % get mean horizontal and vertical gaze position in "clean" interval
            et_hor = mean(EEG_stim_nosac.data([NCHANS_EEG+2 NCHANS_EEG+4],start:stop,t),1);
            et_ver = mean(EEG_stim_nosac.data([NCHANS_EEG+3 NCHANS_EEG+5],start:stop,t),1);
            diff_hor = abs( max(et_hor)-min(et_hor) );
            diff_ver = abs( max(et_ver)-min(et_ver) );
            
            jumpdata   = EEG_stim_nosac.data([NCHANS_EEG+2:NCHANS_EEG+5],start:stop,t);
            jumpdata_d = abs(diff(jumpdata,1,2)); % do first-order diff along 2nd dimension (time)
            maxjump    = max(max(jumpdata_d));
            
            % determine max. gaze change
            if diff_hor > GAZETHRESHOLD_PIXEL || diff_ver > GAZETHRESHOLD_PIXEL || maxjump > MAXGAZEJUMP
                killme = [killme; t];
                warning('killed trial because of gaze deviations');
                counter = counter+1;
            end
        end
        if ~isempty(killme)
            
            %             % plot single-trial gaze (feedback, debugging)
            %             %time = -100:2:4998;
            %             figure;
            %             title(sprintf('dataset: %i subj: %i bad: %i of %i',dataset,s, length(killme), size(EEG_stim_nosac,3) ));
            %             hold on
            %             ylim([200 900]); %xlim([-100 MS_WITHOUT_SAC]);
            %             for j = 1:length(killme)
            %                 %plot(time,EEG_stim_nosac.data(47:50,:,killme(j)));
            %                 plot(EEG_stim_nosac.data(47:50,:,killme(j)));
            %             end
            
            % remove these trials, too.
            EEG_stim_nosac = pop_select(EEG_stim_nosac,'notrial',[killme]);
            EEG_stim_nosac = eeg_checkset(EEG_stim_nosac);
        end
        
        % save single-trial gaze data for both experiments (1 and 2)
        % and for all four ET channels
        if dataset == 1
            singletrialgazeLX1 = [singletrialgazeLX1; squeeze(EEG_stim_nosac.data([NCHANS_EEG+2],start:stop,:))'];
            singletrialgazeLY1 = [singletrialgazeLY1; squeeze(EEG_stim_nosac.data([NCHANS_EEG+3],start:stop,:))'];
            singletrialgazeRX1 = [singletrialgazeRX1; squeeze(EEG_stim_nosac.data([NCHANS_EEG+4],start:stop,:))'];
            singletrialgazeRY1 = [singletrialgazeRY1; squeeze(EEG_stim_nosac.data([NCHANS_EEG+5],start:stop,:))'];
        elseif dataset == 2
            singletrialgazeLX2 = [singletrialgazeLX2; squeeze(EEG_stim_nosac.data([NCHANS_EEG+2],start:stop,:))'];
            singletrialgazeLY2 = [singletrialgazeLY2; squeeze(EEG_stim_nosac.data([NCHANS_EEG+3],start:stop,:))'];
            singletrialgazeRX2 = [singletrialgazeRX2; squeeze(EEG_stim_nosac.data([NCHANS_EEG+4],start:stop,:))'];
            singletrialgazeRY2 = [singletrialgazeRY2; squeeze(EEG_stim_nosac.data([NCHANS_EEG+5],start:stop,:))'];
        end
               
        %% only for EEG_stim_nosac, reject epochs with min/max EEG values (plusminus 120 µV)
        % (can apply strict threshold here because they should contain no ocular artifacts)
        REJECTNOW = 1;
        [EEG_stim_nosac,rejectindex,com] = pop_eegthresh(EEG_stim_nosac, 1, [1:NCHANS_EEG],[-120],[120],EEG.times(1)/1000,MS_WITHOUT_SAC./1000,0,REJECTNOW);
        EEG_stim_nosac = eeg_checkset(EEG_stim_nosac);
                
        warning('\nTrials left with no saccade: %i.\n',size(EEG_stim_nosac.data,3));
        
        %% 2. epochs around saccades
        EEG_sac = pop_epoch(EEG, {'saccade'}, [-0.2 0.6], 'newname', 'saclocked', 'epochinfo', 'yes');
        EEG_sac = applytochannels(EEG_sac,[1:46],'pop_rmbase(EEG,[-100 -10]);');
        EEG_sac = eeg_checkset(EEG_sac);
        
        warning('\nSaccades to analyze found in other epochs: %i.\n',size(EEG_sac.data,3));
        
        %% debugging check: did baseline correction work properly?
        % figure; hold on
        % % plot([-200:2:598],mean(EEG_sac.data(1:45,:,:),3))
        % plot(mean(EEG_sac.data(1:45,:,:),3)')
        % plot(mean(mean(EEG_sac.data([3 4],:,:),1),3)','linewidth',2.0) %IO1+IO2
        % plot(mean(mean(EEG_sac.data([3 4],:,:),1),3)'-mean(mean(EEG_sac.data([33],:,:),1),3)','m','linewidth',1.4) % rEOG
        
        %% get saccades of certain directions (plusminus 30 deg)
        % L-going:    150 - 210 deg
        % R-going:    330 - 360 deg + 0 - 30 deg
        % Up-going:   240 - 300 deg
        % Down-going:  60 - 120 deg
        % Note that the origin of the ET coordinate system is in the upper (!) left corner of the screen...
        
        %% extract saccade angles
        for e = 1:length(EEG_sac.epoch)
            ix = find([EEG_sac.epoch(e).eventlatency{:}] == 0);
            all_sac_angles(e) = cell2mat(EEG_sac.epoch(e).eventsac_angle(ix(1)));
        end
        
        clear ix
        % plot saccade angles
        if PLOT_SAC_ANGLES
            [t,r] = rose(all_sac_angles*pi/180,36); % angle in radians, plot 10° bins
            figure;
            h2 = polar(t,r,'r-');
            set(gca,'ydir','reverse')
        end
        
        %% extract left/right- and up/down-going saccades
        ix_L = find( all_sac_angles > 150  | all_sac_angles < -150); % 150:180 & -150:-180
        ix_R = find( all_sac_angles > -30  & all_sac_angles < 30);   % -30:30
        ix_U = find( all_sac_angles > -120 & all_sac_angles < -60); %  -60:-120 --> because origin of screen is in *upper* left corner!
        ix_D = find( all_sac_angles >   60 & all_sac_angles < 120); %   60:120  --> because origin of screen is in *upper* left corner!

        fprintf('\nNumber of left-right-up-down saccades:')
        [length(ix_L)  length(ix_R) length(ix_U) length(ix_D)]
        clear all_sac_angles % important!
        
        %% make copies of EEG_sac for each direction
        try EEG_sac_L = pop_select(EEG_sac,'trial',ix_L); catch err, EEG_sac_L = []; end
        try EEG_sac_R = pop_select(EEG_sac,'trial',ix_R); catch err, EEG_sac_R = []; end
        try EEG_sac_U = pop_select(EEG_sac,'trial',ix_U); catch err, EEG_sac_U = []; end
        try EEG_sac_D = pop_select(EEG_sac,'trial',ix_D); catch err, EEG_sac_D = []; end
        whos EEG*
        
        %% save benchmark data: EEG_stim, EEG_stim_nosac, EEG_sac_R, // EEG_sac_U, EEG_sac_DEEG_sac_L
        switch dataset
            case 1
                savepath = sprintf('Y:/OPTICA/scenes/benchmarkdata_new/final_epochs/benchmarkdata_subj_%i.mat',s);
            case 2
                savepath = sprintf('Y:/OPTICA/reading/benchmarkdata_new/final_epochs/benchmarkdata_subj_%i.mat',s);
        end
        % save
        save(savepath,'EEG_*','MS_WITHOUT_SAC');
        
        fprintf('\n\n-------------------------------------------------------------');
        fprintf('\nDataset: %i, subject: %i',dataset,s);
        fprintf('\nNumber of epochs in EEG_stim: %i',size(EEG_stim.data,3));
        fprintf('\nNumber of epochs in EEG_sac_R: %i',size(EEG_sac_R.data,3));
        fprintf('\n-------------------------------------------------------------\n\n');
        
        clear EEG*
    end % subject loop
end % dataset loop

%% save single-trial gaze trajectories
save Y:\OPTICA\results\singletrialgaze.mat singletrialgaze* counter
save Y:\OPTICA\results\eyemovements.mat eyemovements

fprintf('\n\nDone.')