% demo_obj.m
% Parallel to demo.m, but driving the whole pipeline through the
% FMRI_Activate_Obj class instead of the low-level fmri_act / fmri_tscore /
% overlayview calls. Uses the same simulated data so the two are comparable.
%
% Requires MIRT on the path (git@github.com:JeffFessler/MIRT.git) 

%% set fMRI parameters
config = struct();
config.t_on  = 20;   % on time (s)
config.t_off = 20;   % off time (s)
config.tr    = 2;    % temporal resolution (s)
config.ncyc  = 5;    % number of cycles to acquire
config.onset_delay = 0;   % paradigm-vs-acquisition offset used by the GLM design (s)
al = 5e-3;   % activation level
nl = 1e-3;   % noise level
N  = 64;     % 2D image matrix size

% time vector
nframes  = floor((config.t_on + config.t_off) * config.ncyc / config.tr);
config.t = config.tr * (0:nframes-1);

%% simulate data (inverse crime): static phantom + boxcar activation + noise
act     = fmri_act(config.t, config.t_on, config.t_off, 0.4); % true response, slight 0.4 s delay
act_msk = phantom([al,.2,.1,0,-0.8,0], N);
und     = phantom(N);
img = und + ...                                 % static component (Shepp-Logan phantom)
      act_msk .* permute(act(:),[2,3,1]) + ...  % dynamic (activation) component
      nl * randn(N, N, nframes);                % noise component

%% build the activation object (constructor builds the GLM design A from config)
fmri_act_obj = FMRI_Activate_Obj(config, img);

%% run the pipeline
% fmri_act_obj.doPolyDetrend();     % optional: remove low-order (drift) trend
fmri_act_obj.get_tscore();          % compute the activation-contrast t-score map
fmri_act_obj.display_tscore();      % show the raw t-score map

%% overlay the activation on the anatomy (first volume as the anatomy layer)
fmri_act_obj.overlay_activation(10, 20);   % t-threshold 10, colormap max 20

% ...or choose a threshold interactively; for 3D volumes the optional 2nd/3rd
% args set the colormap max and the number of montage rows:
% fmri_act_obj.pick_t_thrld(20);        % default 1 row
% fmri_act_obj.pick_t_thrld(20, 2);     % 2 montage rows
