# fMRI-processing
Data analysis for fMRI

`/processing` directory contains code for processing fmri data
- `demo.m` contains a demo fMRI processing pipeline for simulated fMRI data
- `fmri_tscore.m` calculates tscores for a given image timeseries
- `fmri_act.m` creates an hrf-convoled stimulus regressor
- `overlayview.m` displays overlaid activation maps onto a background image using MIRT (git@github.com:JeffFessler/mirt.git)
- `polydetrend.m` performs polynomial detrending on temporal data

`/test-retest-reliability` contains implementation of the statistic model poposed by Genovese et al. for fMRI activation reliability estimation
