function ref = fmri_act(t, t_on, t_off, t_del, varargin)
% function to create fMRI HRF-convolved activation waveform with HRF
% parameters based on defaults from SPM
% by David Frey
%
% inputs:
% t - time vector (s)
% t_on - duration of stimulus "on" period (s)
% t_off - duration of stimulus "off" period (s)
% t_del - stimulus onset delay (s)
%
% HRF parameters (varargin):
% del_resp - response delay (s)
% del_undr - undershoot delay (s)
% dis_resp - response dispersion (s)
% dis_undr - undershoot dispersion (s)
% amp_ratio - response:undershoot amplitude ratio
% onset - total onset delay (s)
%
% outputs:
% ref - reference activation waveform
%

    % set default hrf parameters
    arg.del_resp = 6;
    arg.del_undr = 16;
    arg.dis_resp = 1;
    arg.dis_undr = 1;
    arg.amp_ratio = 6;
    arg.onset = 0;
    arg = vararg_pair(arg, varargin);
    
    % get tr
    tr = t(2) - t(1);
    nt = length(t);
    
    % make hrf
    nk = ceil(32/tr); % kernel length
    resp = gampdf(0:nk, arg.del_resp/arg.dis_resp, arg.dis_resp/tr);
    undr = gampdf(0:nk, arg.del_undr/arg.dis_undr, arg.dis_undr/tr);
    hrf = arg.amp_ratio*resp - undr;
    hrf = hrf/sum(hrf);
    
    % create stimulus waveform
    stim = 1*(mod(t - t_del, t_on + t_off) >= t_off).*(t >= t_del);
    
    % convolve to create activation waveform
    ref = conv(hrf(:).',stim(:).');
    ref = ref(1:nt);

end

