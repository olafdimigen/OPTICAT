%% Simple pipeline to remove ocular artifacts from (free viewing) EEG 
% using optimized ICA training (OPTICAT) including eye tracker-guided 
% component identification (Plöchl et al., 2012) and eye tracker-based
% quality control (under-/overcorrection)
%
%% This script implements procedures from:
% Dimigen, O. (2018). Optimized ICA-based removal of ocular artifacts
% from free viewing EEG. BioArXiv
%
% Please cite this preprint if you use or adapt this script. Thanks!
% olaf.dimigen@hu-berlin.de, Script version: 2018-10-30


%% Load your EEG dataset (continuous or epoched)
% This dataset needs to already include 'saccade' and 'fixation' events in 
% the EEG.event structure. These can be added with the EYE-EEG toolbox 
% (http://www2.hu-berlin.de/eyetracking-eeg). 

%% After loading your data, run the following script:
fprintf('\nCreating optimized ICA training data (OPTICAT)...')

%% Constants
HIPASS           = 2.5  % in Hz. Note: this is the edge of the passband (attenuation of -3 dB). Try even higher value for Reading
OW_FACTOR        = 1    % value for overweighting of SPs
REMOVE_EPOCHMEAN = true % mean-center the appended peri-saccadic epochs? (strongly recommended)
EEG_CHANNELS     = 1:45 % indices of all EEG channels (exclude ET channels here)

%% Create training data & high-pass filter it 
EEG_training = pop_eegfiltnew(EEG,HIPASS,[]); 

%% Cut training data into epochs, e.g. around stimulus onsets (as in Dimigen, 2018)
EEG_training = pop_epoch(EEG_training,{'S123','S234'},[-0.2 2.8]);

%% Remove epoch mean (Groppe, Makeig, & Kutas, 2009)
EEG_training = pop_rmbase(EEG_training,[]);

%% Overweight spike potentials: append peri-saccadic intervals to training data
EEG_training = pop_overweightevents(EEG_training,'saccade',[-0.02 0.01],OW_FACTOR,REMOVE_EPOCHMEAN); % -20 to +10 ms

%% Run ICA
fprintf('\nRunning ICA on optimized training data...')
EEG_training = pop_runica(EEG_training,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % (use BINICA for increased speed)

%% Remember weights & sphering matrix 
wts = EEG_training.icaweights;
sph = EEG_training.icasphere;

%% Remove existing ICA solutions from original dataset
EEG.icaact      = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.icawinv     = [];

%% Transfer unmixing weights
fprintf('\nTransfering ICA weights from training data to original dataset...')
EEG.icasphere   = sph;
EEG.icaweights  = wts;
EEG.icachansind = EEG_CHANNELS;
EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv

fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')

%% Eye-tracker-guided selection of ICs
IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined suitable in Dimigen, 2018)
OPTIMAL_WINDOW   = [5 0]; % saccade window (in samples!) (suitable is -10 ms to 0, see Dimigen, 2018)

PLOTFIG          = true;  % plot figure visualizing influence of threshold?
ICPLOTMODE       = 2;     % plot IC topos (inverse weights)? (2 = only plot ocular ICs)
FLAGMODE         = 3;     % overwrite existing rejection flags?

%% Automatically flag ocular ICs (Plöchl et al., 2012)
[EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade','fixation',OPTIMAL_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);

%% Remove flagged ICs
badcomps = EEG.reject.gcompreject;
EEG      = pop_subcomp(EEG,find(badcomps)); % remove


%% ---------------------------- QUALITY CHECK ----------------------------

%% Evaluate results: Compute saccade-locked ERPs for saccades of a single orientation
EEG_sac = pop_epoch(EEG,{'saccade'},[-0.2 0.6]); % cut epochs around saccades

%% Get saccades of certain orientations (here: plusminus 30 deg)
% (note that the origin of most ET systems is in *upper* left screen corner)
% Right-going:      0 degree
% Left-going:     180 degree
% Up-going:       -90 degree (for most ET data)
% Down-going:      90 degree (for most ET data)

% Extract saccade angles
for e = 1:length(EEG_sac.epoch)
    ix = find([EEG_sac.epoch(e).eventlatency{:}] == 0);
    sac_angles(e) = cell2mat(EEG_sac.epoch(e).eventsac_angle(ix(1)));
end
% plot saccade angles
[t,r] = rose(sac_angles*pi/180,36); % angle in radians, plot 10° bins
figure('name','saccade angles')
h2 = polar(t,r,'r-');
set(gca,'ydir','reverse') % origin of coordinate system upper left

%% Extract left/right/up/down-going saccades (plusminus 30°)
ix_R = find( sac_angles >  -30 & sac_angles <   30);
ix_L = find( sac_angles >  150 | sac_angles < -150);
ix_U = find( sac_angles > -120 & sac_angles <  -60);
ix_D = find( sac_angles >   60 & sac_angles <  120);

%% Cut epochs around saccades (here: rightward)
try 
    fprintf('\nNumber of \"rightwards\" saccades: %i',length(ix_R))
    EEG_sac_R = pop_select(EEG_sac,'trial',ix_R); % select rightwards sacc.
    EEG_sac_R = pop_rmbase(EEG_sac_R,[-100 0])    % remove pre-saccade baseline
catch err
    EEG_sac_R = [];
end

%% Plot saccade-ERP to check for residual artifacts 
% (here: rightward saccades)
figure('name','Control plot: Undercorrection')
plot(EEG_sac_R.times, mean(EEG_sac_R.data(EEG_CHANNELS,:,:),3))
ylabel('Saccade-related ERP')
xlabel('Time after saccade [ms]')

% You can now re-run this script with more traditional settings 
% (e.g. high-pass: 0.5 Hz, no overweighting, low-pass: 40 Hz) and compare 
% the correction results

%% You can also test for overcorrection by ICA (Dimigen, 2018)
% by looking at distortions introduced by ICA in (virtually) 
% eye movement-free intervals without (micro)saccades.
% Code for this will be added later to this script. Or drop me an email
% at olaf.dimigen@hu-berlin.de
