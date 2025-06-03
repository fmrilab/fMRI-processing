%% set fMRI parameters
t_on = 20; % on time (s)
t_off = 20; % off time (s)
tr = 2; % temporal resolution (s)
ncyc = 5; % number of cycles to acquire
al = 5e-3; % activation level
nl = 1e-3; % noise level
N = 64; % 2D image matrix size

% create time vector
nframes = floor((t_on+t_off)*ncyc/tr);
t = tr * (0:nframes-1);

%% simulate data (via inverse crime)
act = fmri_act(t, t_on, t_off, 0.4); % introduce a slight delay
act_msk = phantom([al,.2,.1,0,-0.8,0], N);
und = phantom(N);
img = und + ... % static component (shepp logan phantom)
     act_msk.*permute(act(:),[2,3,1]) + ... % dynamic component
     nl*randn(N,N,nframes); % noise component

%% create the GLM
ref = fmri_act(t, t_on, t_off, 0);
A = ref(:) .^ [0,1]; % model regressors: [baseline, activation]

%% calculate tscores
[tscore, beta] = fmri_tscore(A, abs(img));

%% visualize (requires MIRT: git@github.com:JeffFessler/MIRT.git)
ov = overlayview;
ov.addlayer(und(:,:,1));
ov.addlayer(tscore(:,:,2), ...
    'caxis',[10,20],'cmap','hot','name','tscore');
ov.show
title('Activation map');
axis off