%% create high-pass filtered copies of the data
clear; close all

addpath M:/Dropbox/_subfunc_master % pop_eegfiltnew_ignore_dccorr()
addpath M:/Dropbox/eeglab14_1_1b/
eeglab; close;

% specific constants
LOWCUTOFFS  = [0 0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 5 7.5 10 12.5 15 20 25 30]; % Hz
HICUTOFFS   = [40 100]; % Hz % note: called "LOW/HICUTOFFS" in all scripts, but actually the passband edges
SUBJECTS    = 1:12;

% General note: whenever the scripts say LOWCUTOFFS/HIGHCUTOFFS, this always refers
% to the edge of the passband, since this is what pop_eegfiltnew() expects as input
% (rather than the [half-amplitude/-6dB] cutoff value)

% dataset loop
for dataset = 1
    
    switch dataset
        case 1 % scenes
            path_rawset_et  = 'Y:/OPTICA/scenes/raw_set_et/';
            path_filt       = 'Y:/OPTICA/scenes/filt_new/';
        case 2 % reading
            path_rawset_et  = 'Y:/OPTICA/reading/raw_set_et/';
            path_filt       = 'Y:/OPTICA/reading/filt_new/';
    end
    
    % subject loop
    for s = SUBJECTS
        subtime = tic;
        
        fprintf('\n\n\n---------------------------------------------');
        fprintf('\nProcessing dataset %i, subject %i',dataset,s)
        fprintf('\n---------------------------------------------\n');
                        
        % low-pass filter loop
        for hicutoff = HICUTOFFS
            
            clear EEG*
            
            %% load original recording (with zero-A1 channel and synchronized ET)
            loadfilename = sprintf('eeg_%02i.set',s);
            EEG = pop_loadset('filename',loadfilename,'filepath',path_rawset_et);
                       
            %% remove the mastoid online reference (A1) and ET data (reappended later)
            EEG_tmp      = EEG.data(46:end,:);
            chanlocs_tmp = EEG.chanlocs(46:end);
            
            EEG.chanlocs(46:end) = [];
            EEG.data(46:end,:)   = [];
            EEG.nbchan = 45;
            
            % EEG = eeg_checkset(EEG);
            % eeglab redraw
                                    
            if hicutoff < 100 % don't filter at 100 Hz // already done in raw data
                fprintf('\n\n######################################');
                fprintf('\nFiltering with High Cutoff of: %.1f Hz...',hicutoff);
                fprintf('\n######################################\n\n');
                EEG = pop_eegfiltnew(EEG,[],hicutoff);
            end
            
            % create copy that will be high-pass filtered differently
            EEG2 = EEG; clear EEG
            
            % high-pass filter loop
            for lowcutoff = LOWCUTOFFS
                
                if lowcutoff == 0 % do nothing, raw data already filtered at 10 sec TC
                    fprintf('\n\n######################################');
                    fprintf('\nNot filtering, keeping 0.016 Hz filter')
                    fprintf('\n######################################\n\n');
                    EEG = EEG2;
                else
                    fprintf('\n\n######################################');
                    fprintf('\nFiltering with Low Cutoff of %.1f Hz...',lowcutoff);
                    fprintf('\n######################################\n\n');
                    
                    switch dataset
                        case 1
                            EEG = pop_eegfiltnew(EEG2,lowcutoff,[]); % start filtering 400 samples after DC corr
                        case 2
                            % reading dataset was recorded as DC data
                            % this means that the amp sometimes does a DC correction that resets all channels
                            % creating an artifact that can distort a LP-filter. We use a slightly
                            % modified version of pop_eegfiltnew that
                            % separately filters the intervals in-between possible DC_corrections
                            % filter between data breaks (if any)
                            EEG = pop_eegfiltnew_ignore_dccorr(EEG2,lowcutoff,[],400); % start filtering 400 samples after DC corr
                    end
                end
                
                % re-append ET data and channels A1
                EEG.data(46:50,:)   = EEG_tmp;
                EEG.chanlocs(46:50) = chanlocs_tmp;
                EEG.nbchan = 50;
                EEG.setname = ['dataset_' num2str(dataset) '_sub_' num2str(s) '_freq_' num2str(lowcutoff) '_' num2str(hicutoff)];
                EEG = eeg_checkset(EEG);
                % eeglab redraw
                
                % save filtered dataset
                savefilename = sprintf('eeg_%02i_lowcut_%.1fHz_hicut_%iHz.set',s,lowcutoff,hicutoff);
                EEG = pop_saveset(EEG,'filename',savefilename,'filepath',path_filt);
                
                clear EEG
                
            end % highpass loop
        
        end % lowpass loop
        
        t_for_subject = toc(subtime);
        warning('Time info following:')
        fprintf('\n\n-----------------------------------------------');
        fprintf('\nSeconds for subject %i, dataset: %i: %.2f',s,dataset,t_for_subject);
        fprintf('\n-----------------------------------------------\n\n');
        
    end % subject
end % datasets

fprintf('\n\nDone')