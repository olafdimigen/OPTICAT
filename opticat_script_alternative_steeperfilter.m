% Alternative version of the "OPTICAT" script that uses a "steeper" 
% high-pass filter with a custom, more narrow transition bandwidth
%
% Explanation: Procedures in the Dimigen et al. (2020) paper use 
% EEGLAB's pop_eegfiltnew() function with its default transition bandwidth 
% settings. Per default, this function filters the EEG using a relatively wide
% and frequency-dependent transition bandwidth 
% (see https://github.com/downloads/widmann/firfilt/firfilt.pdf).
%
% Removing low frequencies from the training data likely means that
% low-frequency neural activity (e.g. < 1-2 Hz) will not be properly modelled 
% in the ICA decomposition. In the paper, I did not find that this had a
% problematic effect on the stimulus-ERPs, that is, I did not see any 
% increase in overcorrection of the stimulus-ERP if the *training
% data* was high pass-filtered aggressively (e.g. at 2-2.5 Hz) with the 
% EEGLAB defaults.

% However, Supplementary Figure S7 suggests that the problematic frequency
% content in the training data is below 1 Hz, at least for Scene viewing.
% Thus, a potentially safer alternative is to filter the EEG at a lower cutoff
% (e.g. at 0.75 Hz, see Supplementary Figure S7) but with a "steeper" filter.
% 
% This should remove less of the low-frequency neural activity
% from the training data while producing comparable correction results 
% (cf. Supplementary Figure S7) as the optimal filter settings described with 
% the EEGLAB filter defaults detailed in the paper. Preserving low-frequency 
% activity in the training data would be beneficial if the user wants to work 
% with the decomposed source waveforms in source space (rather than backproject 
% the corrected data).
%
% The script below implements the OPTICAT procedures with a custom, 
% "steeper" filter setting. For more documentation on the OPTICAT procedure
% itself, see the "opticat_script.m" in this folder.
%
% I would like to thank an anonymous contributor for helpful input
% on filter design.

clear

%% add the firfilt toolbox by A. Widmann (https://github.com/widmann/firfilt)
addpath M:\Dropbox\tools\firfilt-master

%% Load your EEG dataset (can be continuous or epoched)
% Note: This dataset in EEGLAB format needs to already include 
% 'saccade' and 'fixation' events in the EEG.event structure.
% These can be added, for example, with the EYE-EEG toolbox. 
EEG = pop_loadset('filename','C:/myEEGdata.set');


OW_FACTOR        = 1    % value for overweighting of SPs (1 = add data corresponding to 100% of original data length)
REMOVE_EPOCHMEAN = true % mean-center the appended peri-saccadic epochs? (strongly recommended)
EEG_CHANNELS     = 1:45 % indices of all EEG channels (exclude any eye-tracking channels here)
                        % I recommend to also include EOG channels if they
                        % were also recorded against the common reference % (not: bipolar)

% ------------------START: STEEPER HIGH-PASS FILTER -----------------------
CUTOFF           = 0.75 % -6 dB cutoff (in Hz)
                        % Note: possibly go even lower, e.g. to 0.5 Hz if your priority is
                        % to preserve low-frequency neural activity in the decomposition. 
                        % However, compare Suppl. Figure S7 for effects on correction quality
TBW              = 0.5  % transition bandwidth (in Hz)

% compute filter order
filterorder = firwsord('hamming', EEG.srate, TBW);

% design windowed sinc filter
b = firws(filterorder, CUTOFF / (EEG.srate/2), 'high', windows('hamming', filterorder+1));

% get frequency response
[h,f] = freqz(b,1,16384,EEG.srate);

% visualize filter characteristics
figure
subplot(2,1,1); hold on; title('linear')
plot(f,abs(h))
xlim([0 2]), ylim([0 1.05])
subplot(2,1,2); hold on; title('dB')
plot(f,20*log10(abs(h)))
xlim([0 2]), ylim([-80 0])
xlabel('Hz')
% -------------------END: STEEPER HIGH-PASS FILTER ------------------------


%% Create copy of data used as training data & high pass-filter it
EEG_training = EEG;
EEG_training = firfilt(EEG_training,b,[]); % filter the training data using firfilt()

%% Cut training data into epochs, e.g. around stimulus onsets (as in Dimigen, 2020)
EEG_training = pop_epoch(EEG_training,{'S 11','S 12','S 13'},[-0.2 1.8]);

%% Overweight spike potentials
% Repeatedly append intervals around saccade onsets (-20 to +10 ms) to training data
EEG_training = pop_overweightevents(EEG_training,'saccade',[-0.02 0.01],OW_FACTOR,REMOVE_EPOCHMEAN);

%% Run ICA
fprintf('\nRunning ICA on optimized training data...')
EEG_training = pop_runica(EEG_training,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % or use binary ICA for more speed

%% Remember ICA weights & sphering matrix
wts = EEG_training.icaweights;
sph = EEG_training.icasphere;

%% Remove any existing ICA solutions from your original dataset
EEG.icaact      = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.icawinv     = [];

%% Transfer unmixing weights
fprintf('\nTransfering ICA weights from training data to original data...')
EEG.icasphere   = sph;
EEG.icaweights  = wts;
EEG.icachansind = EEG_CHANNELS;
EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv


fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')

%% Eye-tracker-guided selection of ICs
IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined as suitable in Dimigen, 2020)
SACC_WINDOW      = [5 0]; % saccade window (in samples!) to compute variance ratios (see Dimigen, 2020)
PLOTFIG          = true;  % plot a figure visualizing influence of threshold setting?
ICPLOTMODE       = 2;     % plot component topographies (inverse weights)? (2 = only plot "bad" ocular ICs)
FLAGMODE         = 3;     % overwrite existing rejection flags? (3 = yes)

%% Automatically flag ocular ICs (Plöchl et al., 2012)
[EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade','fixation',SACC_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);

%% Remove flagged ocular ICs
badcomps = EEG.reject.gcompreject;
EEG      = pop_subcomp(EEG,find(badcomps)); % remove them



%%  ----- SIMPLE CORRECTION QUALITY CHECK (TEST FOR UNDERCORRECTION) -----

%% Evaluate results: Compute saccade-locked ERPs for saccades of a single orientation
EEG_sac = pop_epoch(EEG,{'saccade'},[-0.2 0.6]); % cut epochs around saccades

%% Get saccades of certain orientations (here: plusminus 30 deg)
% (note that the origin of most eye-tracker systems is *upper* left screen corner)
% Right-going:      0 degree
% Left-going:     180 degree
% Up-going:       -90 degree (for most ET systems)
% Down-going:      90 degree (for most ET systems)

% Extract saccade angles
for e = 1:length(EEG_sac.epoch)
    ix = find([EEG_sac.epoch(e).eventlatency{:}] == 0);
    sac_angles(e) = cell2mat(EEG_sac.epoch(e).eventsac_angle(ix(1)));
end

% Plot saccade angles
[t,r] = rose(sac_angles*pi/180,36); % angle in radians, plot 10° bins
figure('name','saccade angles')
h2 = polar(t,r,'r-');
set(gca,'ydir','reverse') % origin of coordinate system = upper left corner

%% Extract left/right/up/down-going saccades (plusminus 30°)
ix_R = find( sac_angles >  -30 & sac_angles <   30);
ix_L = find( sac_angles >  150 | sac_angles < -150);
ix_U = find( sac_angles > -120 & sac_angles <  -60);
ix_D = find( sac_angles >   60 & sac_angles <  120);

%% Cut epochs around saccades (here: rightward sacc. only, adapt as you wish)
try
    fprintf('\nNumber of \"rightwards\" saccades: %i',length(ix_R))
    EEG_sac_R = pop_select(EEG_sac,'trial',ix_R); % select rightwards sacc.
    EEG_sac_R = pop_rmbase(EEG_sac_R,[-100 0],[]);   % subtract pre-sacc. baseline
catch err
    EEG_sac_R = [];
end

%% Plot saccade-ERP to check for residual artifacts
% (here: for rightward saccades)
figure('name','Control plot: Undercorrection')
plot(EEG_sac_R.times, mean(EEG_sac_R.data(EEG_CHANNELS,:,:),3))
ylabel('Saccade-related ERP')
xlabel('Time after saccade [ms]')
