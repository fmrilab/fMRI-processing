function tc = mkref2(time, ontimes, offtimes)
% times in s, tr = time(2) - time(1)
% vector of ontimes and off times in s (same length)


% nslices*scanner_tr and convert to s
tr = time(2) - time(1); % temporal resolution of fMRI data (s)
ntime = length(time); % # of time points 
totaltime = ntime*tr;

% make hrf
hrf = spm_hrf(tr);

% define task in hrftr's
task = zeros(size(time));
for lp = 1:length(ontimes)
    task = task + (time >= ontimes(lp)).*(time < offtimes(lp)); % task times in s
end

ref1 = conv(hrf,task);

% resample to tr
tvec2 = [1:length(ref1)]*tr;

tc = ref1(1:length(time));