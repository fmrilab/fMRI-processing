# fMRI-processing
Data analysis for fMRI

`/processing` directory contains code for processing fmri data
- `fmri_tscore.m` calculates tscores for a given image timeseries
- `mkref2.m` creates an hrf-convoled stimulus regressor

`/visualization` directory contains code for visualizing data
- `overlayview.m` displays overlaid activation maps onto a background image

`/test-retest-reliability` contains implementation of the statistic model poposed by Genovese et al. for fMRI activation reliability estimation
