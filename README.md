Create optimized ICA training (OPTICAT) data for ocular correction of EEG recorded during free viewing (or data that is otherwise contaminated by ocular artifacts)

This repository contains two things:

1. A simple script to run an optimized eye tracker-guided ICA on your own data
2. Code of the main analyses in the corresponding paper

Reference paper: 

Dimigen, O. (2020). Optimizing the ICA-based removal of ocular artifacts
from free viewing EEG. 116117. NeuroImage, https://doi.org/10.1016/j.neuroimage.2019.116117 


Please see the scripts itself for more infos und documentation.

Olaf Dimigen
olaf.dimigen@hu-berlin.de, Jan. 2020


Minor update 2024-06-27:

A few users asked whether the OPTICAT method also flags the eye-lid artifacts caused by blinks. This is often the case if the training data contains enough vertical eye movements,
since vertical saccades also go along with a small blink-like activity (see Introduction section of the paper). However, in other experiments, e.g. those only involving 
only horizontal saccades, it is possible that classic blink components are not flagged by the variance-ratio criterion. 

If you encounter problems with OPTICAT not modeling or not flagging blink components, I suggest the following steps:

1. Create optimized ICA training data as described in the paper (but make sure the training data also contains some blinks!)
2. Run ICA on it
3. Flag the "bad" components with OPTICAT and "remember" the component indices of the bad ICs
4. Additionally use EEGLAB's "ICLABEL" plugin on this optimized ICA decomposition, with only the "EYE" criterion activated
(note: while ICLABEL is not always good at flagging ocular artifacts like the spike potential, it will reliably flag blinks)
5. Combined the bad component indices from both 3. and 4. and remove these ICs

