%% Pipeline to remove ocular artifacts from (free viewing) EEG 
% using optimized ICA training (OPTICAT, Dimigen, 2018) and 
% eye tracker-guided component identification (Plöchl et al., 2012)
%
%% This pipeline implements procedures from:
% Dimigen, O. (2018). Optimizing ICA-based ocular artifact correction 
% of free viewing EEG. BioArXiv
%
% Please cite this preprint if you use or adapt this script. Thanks!
% olaf.dimigen@hu-berlin.de 
% Script version: 2018-10-22

%% Load your EEG dataset (continuous or epoched)
% This dataset needs to already include 'saccade' and 'fixation' events in 
% EEG.event. These can be added with the EYE-EEG toolbox 
% (http://www2.hu-berlin.de/eyetracking-eeg). 
%
% (It's probably also possible to run the procedures with saccades detected
% in a radial rEOG channel)
% EEG = pop_loadset('filename','my_eegfilename.set','filepath','my_filepath');

%% After loading the data, run the following script:
fprintf('\nCreating optimized ICA training data (OPTICAT)...')

%% Constants
HIPASS           = 2.5  % Hz
OW_PROPORTION    = 1    % value for overweighting of SPs used in Dimigen (2018)
REMOVE_EPOCHMEAN = true % mean-center the appended peri-saccadic epochs? (strongly recommended)
EEG_CHANNELS     = 1:45 % indices of EEG channels (exclude ET channels here)

%% Create training data and high-pass filter it 
EEG_opticat = pop_eegfiltnew(EEG,HIPASS,[]); 

%% Cut training data into epochs, e.g. around stimulus onsets (as in Dimigen, 2018)
% EEG_opticat = pop_epoch(EEG_opticat,{'S123','S234'},[-0.2 2.8]);

%% Remove epoch mean for existings epochs (Groppe, Makeig, & Kutas, 2009)
EEG_opticat = pop_rmbase(EEG,[]);

%% Overweight peri-saccadic intervals (containing spike potentials)
EEG_opticat = pop_overweightevents(EEG_opticat,'saccade',[-0.02 0.01],OW_PROPORTION,REMOVE_EPOCHMEAN);

%% Run ICA
fprintf('\nRunning ICA on optimized training data...')
EEG_opticat = pop_runica(EEG_opticat,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % or use BINICA for increased speed

%% Remember weights & sphering matrix 
wts = EEG_opticat.icaweights;
sph = EEG_opticat.icasphere;

%% Remove existing ICA solutions from original dataset
EEG.icaact      = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.icawinv     = [];

%% Transfer "optimized" weights
fprintf('\nTransfering ICA weights from training data to original dataset...')
EEG.icasphere   = sph;
EEG.icaweights  = wts;
EEG.icachansind = EEG_CHANNELS;
EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv

fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')

%% Eye-tracker-guided selection of ICs
IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined suitable in Dimigen, 2018)
OPTIMAL_WINDOW   = [5 0]; % saccade window (in samples) (determined as suitable in Dimigen, 2018)
PLOTFIG          = true;  % plot figure visualizing influences of the IC_THRESHOLD?
ICPLOTMODE       = 2;     % plot IC topos (inverse weights)? Here: only plot "bad" ocular ICs
FLAGMODE         = 3;     % overwrite existing flags

%% Automatically flag ocular ICs (Plöchl et al., 2012)
[EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade','fixation',OPTIMAL_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);

%% Remove flagged ICs
badcomps = EEG.reject.gcompreject;
EEG = pop_subcomp(EEG,find(badcomps));

%% Evaluate results: Compute saccade-locked ERPs for saccades of single orientation
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
    EEG_sac_R = pop_rmbase(EEG_sac_R,[-200 0])    % remove pre-saccade baseline
catch err
    EEG_sac_R = []
end

%% Plot saccade-ERP to check for residual artifacts 
% here: only rightward saccades
figure('name','Ocular correction control plot')
plot(EEG_sac_R.times, mean(EEG_sac_R.data(EEG_CHANNELS,:,:),3))
ylabel('Saccade-related ERP')
xlabel('Time after saccade [ms]')

% for comparison, re-run this script with more traditional settings 
% (e.g. high-pass: 0.5 Hz, no overweighting, low-pass: 40 Hz)

%% You can also test for overcorrection by ICA (Dimigen, 2018)
% by looking at distortions introduced by ICA in (virtually) 
% eye movement-free intervals without (micro)saccades.
