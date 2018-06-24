%% Optimized pipeline to remove ocular artifacts from (free viewing) EEG 
% including Optimized ICA Training data (OPTICAT, Dimigen, 2018) and 
% eye tracker-guided identification of ocular ICs (Plöchl et al., 2012)
%
%% This pipeline is based on the following preprint:
%
% Dimigen, O. (2018). Optimizing ICA-based ocular artifact correction 
% of EEG data recorded during free viewing. BioArXiv
%
% Please cite this preprint if you use or adapt this script. Thanks!
%
% olaf.dimigen@hu-berlin.de, June 2018


%% Load your continuous or epoched EEG dataset
% this dataset should already have 'saccade' and 'fixation' events 
% detected with EYE-EEG in the EEG.event structure. 
%
% EEG = pop_loadset('filename','your_eegfilename.set','filepath','your_filepath');


%% After loading the data, run the following script:
fprintf('\nCreating optimized ICA training data (OPTICAT)...')

%% Constants
HIPASS           = 2.5   % Hz
OW_PROPORTION    = 1     % value for overweighting of SPs used in Dimigen (2018)
REMOVE_EPOCHMEAN = true  % subtract mean from appended epochs? (recommended)
EEG_CHANNELS     = 1:45  % indices of EEG channels (exclude ET channels etc.)

%% Make training data and high-pass filter it 
EEG_opticat = pop_eegfiltnew(EEG,HIPASS,[]); 

%% Remove epoch mean for existings epochs (Groppe, Makeig, & Kutas, 2009)
EEG_opticat = pop_rmbase(EEG,[]);

%% Overweight intervals around saccades (containing spike potential)
EEG_opticat = pop_overweightevents(EEG_opticat,'saccade',[-0.02 0.01],OW_PROPORTION,REMOVE_EPOCHMEAN);

%% Run ICA
fprintf('\nRunning ICA on the optimized data...')
EEG_opticat = pop_runica(EEG_opticat,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % or use BINICA for increased speed

%% Remember weights & sphering matrix from optimized dataset
wts = EEG_opticat.icaweights;
sph = EEG_opticat.icasphere;

%% Remove existing ICA solutions from original dataset
EEG.icaact      = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.icawinv     = [];

%% Transfer "optimized" solution
fprintf('\nTransfering ICA weights for optimized training data to original dataset...')
EEG.icasphere   = sph;
EEG.icaweights  = wts;
EEG.icachansind = EEG_CHANNELS;
EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv

fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')

%% Eye-tracker-guided selection of ICs
IC_THRESHOLD     = 1.1;   % variance ratio threshold, determined as near-optimal in Dimigen, 2018
OPTIMAL_WINDOW   = [5 0]; % saccade window (in samples), determined as near-optimal in Dimigen, 2018
PLOTFIG          = true;  % plot figure visualizing the IC_THRESHOLD
ICPLOTMODE       = 2;     % plot only topos (inverse weights) of "bad" flagged ICs
FLAGMODE         = 3;     % overwrite existing flags


%% Automatically flag ocular ICs (Plöchl et al., 2012)
[EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade','fixation',OPTIMAL_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);

%% Remove flagged ICs
badcomps = EEG.reject.gcompreject;
EEG = pop_subcomp(EEG,find(badcomps));


%% Optional: Evaluate results: Compute saccade-locked ERP for saccades of single orientation
EEG_sac = pop_epoch(EEG,{'saccade'},[-0.2 0.5]);

% Get saccades of certain orientations (here: plusminus 30 deg)
% (note that origin of most ET systems is in *upper* left screen corner)
% L-going:    150 - 210 deg
% R-going:    330 - 360 deg and 0 - 30 deg
% Up-going:   240 - 300 deg
% Down-going:  60 - 120 deg
        
%% Extract saccade angles
for e = 1:length(EEG_sac.epoch)
    ix = find([EEG_sac.epoch(e).eventlatency{:}] == 0);
    sac_angles(e) = cell2mat(EEG_sac.epoch(e).eventsac_angle(ix(1)));
end
% plot saccade angles
[t,r] = rose(sac_angles*pi/180,36); % angle in radians, plot 10° bins
figure('name','saccade angles')
h2 = polar(t,r,'r-');
set(gca,'ydir','reverse') % origin of coordinate system upper left

%% Extract left/right/up/down-going saccades
ix_R = find( sac_angles >  -30 & sac_angles <   30);
ix_L = find( sac_angles >  150 | sac_angles < -150);
ix_U = find( sac_angles > -120 & sac_angles <  -60);
ix_D = find( sac_angles >   60 & sac_angles <  120);
fprintf('\nNumber of \"rightwards\" saccades: %i',length(ix_R))

%% Cut epochs around rightward saccades
try 
    EEG_sac_R = pop_select(EEG_sac,'trial',ix_R); % select rightwards sacc.
    EEG_sac_R = pop_rmbase(EEG_sac_R,[-200 0])    % remove pre-saccade baseline
catch err
    EEG_sac_R = []
end

%% Plot saccade-ERP: check for residual artifacts
figure('name','Ocular correction control plot')
plot(EEG_sac_R.times, mean(EEG_sac_R.data(EEG_CHANNELS,:,:),3))
ylabel('Saccade-related ERP')
xlabel('Time after saccade')

%% You can also test for overcorrection of the data by ICA (Dimigen, 2018)
% by looking at distortions introduced by ICA in (virtually) 
% eye movement-free intervals