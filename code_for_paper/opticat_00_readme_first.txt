This is a (fairly) cleaned-up version of a series of scripts implementing the main
analyses presented in the paper:

Dimigen (2020). "Optimizing the ICA-based removal of ocular EEG artifacts from free viewing experiments" 
as published in "NeuroImage", 116117 
https://doi.org/10.1016/j.neuroimage.2019.116117 (open access)

All code is in MATLAB except for the (ANOVA) statistics which were run in "R" using the "ez" package.

The scripts are a simplified/cleaned-up version of the pipeline used for the main analyses. 
To keep the pipeline readable, scripts for many of the supplementary analyses (more overweighting proportions, other saccade directions, 
visualization of filter characteristics, control analysis with steeper filter etc.) as well as scripts for some of the 
more fancy visualizations are not included here. Drop me an email if you require those or subfunctions not included here.
The scripts make use of EEGLAB (including the binary ICA implementation) and the EYE-EEG toolbox (v0.85). 
Furthermore, a few functions of the "unfold" toolbox (www.unfoldtoolbox.org) are used for (less important) plotting purposes. 
The (very simple) rmANOVA stats conducted in R require the "EZ" package.

Aggregated subject-level data for all subjects (e.g. ICA inverse weights of all ICA solutions and quality critera) is found at the corresponding OSF repository (https://osf.io/ka6qh/).
I hope the code is helpful for researchers planning similar comparisons. I cannot share the raw EEG/ET data of both experiments at the moment but hope this will be possible in the future.

If you spot any issues or if you have comments or questions, you can contact me at olaf.dimigen@hu-berlin.de

Olaf Dimigen
olaf.dimigen@hu-berlin.de 
Berlin, January 2020 (updated)
