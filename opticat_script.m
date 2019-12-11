%% Pipeline to remove ocular artifacts from (free viewing) EEG 
% using optimized ICA training (OPTICAT), automatic eye tracker-guided 
% component identification, and eye tracker-based quality control
%
% This script implements the procedures from:
% Dimigen, O. (2019). Optimizing the ICA-based removal of ocular artifacts
% from free viewing EEG. NeuroImage, https://doi.org/10.1016/j.neuroimage.2019.116117 
%
% Please cite this publication if you use/adapt this script. Thanks!
%
% For this script to work, you need to have the EYE-EEG extension for
% EEGLAB installed (http://www2.hu-berlin.de/eyetracking-eeg)
% It is also available at: www.github.com/olafdimigen/eye-eeg
%
% olaf.dimigen@hu-berlin.de, Script version: 2019-11-26

%% Constants
HIPASS           = 2    % Filter's passband edge (in Hz)
                        % Best results for scenes were obtained with values of 2 to 2.5 Hz
                        % Possibly try even higher value for tasks like Reading
OW_FACTOR        = 1    % value for overweighting of SPs (1 = add spike potentials corresponding to 100% of original data length)
REMOVE_EPOCHMEAN = true % mean-center the appended peri-saccadic epochs? (strongly recommended)
EEG_CHANNELS     = 1:45 % indices of all EEG channels (exclude any eye-tracking channels here)
                        % I recommend to also include EOG channels (if also recorded against common reference)

%% Load your EEG dataset (can be continuous or epoched)
% Note: This dataset in EEGLAB format needs to already include 
% 'saccade' and 'fixation' events in the EEG.event structure.
EEG = pop_loadset('filename','C:/myEEGdata.set');
     
%% Create copy of data used as training data & high pass-filter it
fprintf('\nCreating optimized ICA training data...')
EEG_train = EEG;
EEG_training = pop_eegfiltnew(EEG_train,HIPASS,[]); 

%% Cut training data into epochs, e.g. around stimulus onsets (as in Dimigen, 2019)
% here I assume stimulus onset triggers are "S123" and "S234"
EEG_training = pop_epoch(EEG_training,{'S123','S234'},[-0.2 2.8]);

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
IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined as suitable in Dimigen, 2019)
SACC_WINDOW      = [5 0]; % saccade window (in samples!) to compute variance ratios (see Dimigen, 2019)
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
    EEG_sac_R = pop_rmbase(EEG_sac_R,[-100 0]);   % subtract pre-sacc. baseline
catch err
    EEG_sac_R = [];
end

%% Plot saccade-ERP to check for residual artifacts 
% (here: for rightward saccades)
figure('name','Control plot: Undercorrection')
plot(EEG_sac_R.times, mean(EEG_sac_R.data(EEG_CHANNELS,:,:),3))
ylabel('Saccade-related ERP')
xlabel('Time after saccade [ms]')

% You can now re-run this script with more traditional settings 
% (e.g. high-pass: 0.5 Hz, low-pass: 40 Hz, no overweighting) and compare 
% the correction results

%% You can also test for overcorrection by ICA (Dimigen, 2019)
% by looking at distortions introduced by ICA in (virtually) 
% eye movement-free epochs without (micro)saccades or large drift movements.
% (I might add this to the script later)

% If this scripts works very nicely or very badly for you, drop me an email
% at olaf.dimigen@hu-berlin.de