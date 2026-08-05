classdef FMRI_Activate_Obj < handle

    properties
        config
        data
    end

    methods
        function obj = FMRI_Activate_Obj(config,img)
            % inputs:
            %   config:
            %       - t_on: on time (s)
            %       - t_off: off time (s)
            %       - tr: temporal resolution (s)
            %       - ncyc: number of cycles
            %       - onset_delay: offset between the paradigm begin and acquisition begin (s)
            %   img: [nx,ny,(nz),nt], the image series (2D or 3D volumes by time)
            obj.config = config;
            obj.initializeData(img);
        end

        function initializeData(obj,img)
            % data:
            %   - img: [nx,ny,(nz),nt]
            %   - A: [nT,2]; model regressor,
            %       A(:,1) : basedline;
            %       A(:,2) : activation;
            %   - tmap: [nx,ny,(nz)] tscore map (activation contrast)
            ref = fmri_act(obj.config.t, obj.config.t_on, obj.config.t_off, obj.config.onset_delay);
            A = ref' .^ [0,1]; %model regressors: [baseline, activation]; the first few are dummy cycles
            obj.data = struct('img',img, ...
                'A',A, ...
                'tmap',[]);
        end

        function doPolyDetrend(obj)
            % detrend DC term due to drifting
            obj.data.img = poly_detrend(obj.data.img,1);
        end

        function get_tscore(obj)
            [tmp,~] = fmri_tscore(obj.data.A,abs(obj.data.img));
            % tmp is [spatial..., ncontrast]; take the activation contrast (2nd
            % regressor) keeping all spatial dims -> works for 2D or 3D volumes.
            sp = repmat({':'}, 1, ndims(tmp)-1);
            obj.data.tmap = tmp(sp{:}, 2);
        end
        
        function display_tscore(obj)
            figure;
            im(obj.data.tmap,'tscore','cbar')
        end

        function overlay_activation(obj,t_thrld,t_max,nrow)
            % overlay t-score map on top of the anatomy
            %   t_thrld: scalar, t-score threshold
            %   t_max(optional): colormap maximum for the tscore map
            %   nrow(optional): number of rows in the overlay montage (default 1)
            if nargin<3 || isempty(t_max), t_max = max(obj.data.tmap(:)); end
            if nargin<4 || isempty(nrow), nrow = 1; end

            viewargs={'row',nrow,'notick'};
            ov = overlayview;
            % first volume (t=1) as the anatomy layer; keep all spatial dims (2D or 3D)
            nsp = ndims(obj.data.img) - 1;
            sp  = repmat({':'}, 1, nsp);
            img = abs(squeeze(obj.data.img(sp{:}, 1)));
            % mask out activation out of the brain
            mask = abs(img)>0.05*max(img(:));

            ov.addlayer(img,'caxis',[min(img(:)),max(img(:))])
            ov.addlayer(obj.data.tmap.*mask,'caxis',[t_thrld,t_max],'cmap','hot');

            figure;
            ov.show('viewargs',viewargs);title('Activation map','FontSize',13)

            obj.config.t_thrld = t_thrld;
        end

        function pick_t_thrld(obj,t_max,nrow)
            % pick a good t-score threshold for overlay view
            % (optional) t_max: colormap maximum for the tscore map display
            % (optional) nrow: number of rows in the overlay montage (default 1)
            % fmri_act_obj.pick_t_thrld(16)        % 1 row 
            % fmri_act_obj.pick_t_thrld(16, 4)     % t_max=16, 4 rows
            % fmri_act_obj.pick_t_thrld([], 4)     % default t_max, 4 rows

            if nargin<2, t_max = []; end
            if nargin<3, nrow = []; end
            t_thrld = input('Give an initial t-score threshold:','s');
            while ~isempty(t_thrld)
                obj.overlay_activation(str2double(t_thrld),t_max,nrow);
                t_thrld = input('Pick a new t-score threshold (hit enter if the value now is good):','s');
            end
        end

    end

end