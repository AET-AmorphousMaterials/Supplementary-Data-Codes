classdef RESIRE_reconstructor
    
    properties
        % input data
        InputProjections
        InputAngles
        Support
        Mask
        InitialModel
        
        % reconstruction and parameters
        reconstruction
        NumProjs
        Dim1  % projection dimension 1
        Dim2  % projection dimension 2
        obj_dimx % object dimension 1
        obj_dimy % object dimension 2
        obj_dimz = 1; % object thickness, 1 - same as projectoin dim1; 2 - same as projection dim2; Other - specified size
        
        n1_oversampled  % projection dimension 1 after oversampling
        n2_oversampled  % projection dimension 2 after oversampling
        n3_oversampled  % projection dimension 3 after oversampling
        
        errR
        Rarr_record
        Rarr2_record
        
        % R factor monitoring
        monitor_R = 1
        monitorR_loopLength = 1
        
        % filenames
        filename_Projections = ''
        filename_Angles = '';
        filename_Support = '';
        filename_Mask = '';
        filename_InitialModel = '';
        filename_Results = './results/RESIRE_rec.mat';
        
        % reconstruction parameters
        numIterations = 50;
        oversamplingRatio_x = 3;
        oversamplingRatio_y = 3;
        oversamplingRatio_z = 3;
        griddingMethod = 1;
        extenedobject = 'false';
        positivity = 'true';
        method = 'FST'; % forward projection method: 'FST' or 'Radon'
        
        % axes vectors for phi, theta, psi rotations
        vector1 = [0 0 1];
        vector2 = [0 1 0];
        vector3 = [1 0 0];
        
        step_size = 1;         % normalized step_size is chosen in range (1,3)
        sigma  = 0;           % smoothing constraint, smoother with smaller sigma, 0 if no smooth
        dtype  = 'single';      % data type, single for memory efficiency
        kernel;
        Plot_rec = 0;
        
        sum_rot_pjs;
    end
    
    methods
        
        % declare long methods in external files
        obj = reconstruct(obj) % FST method
        
        obj = reconstruct3D2D(obj) % Radon transform for single y-tilt case
        
        obj = reconstruct_control(obj)
        % declare short methods in this file
        function obj = readFiles(obj)
            if FileExist(obj.filename_Projections)
                obj.InputProjections = importdata(obj.filename_Projections);
            else
                error('RESIRE: Projections file does not exist!')
            end
            if FileExist(obj.filename_Angles)
                obj.InputAngles = importdata(obj.filename_Angles);
            else
                error('RESIRE: Angles file does not exist!')
            end
            if isempty(obj.filename_Support)
                fprintf('RESIRE: No support used.\n');
            else
                if FileExist(obj.filename_Support)
                    obj.Support = importdata(obj.filename_Support);
                    fprintf('RESIRE: Support used.\n');
                else
                    error('RESIRE: Support file does not exist!')
                end
            end
            if isempty(obj.filename_Mask)
                fprintf('RESIRE: No Mask used.\n');
                obj.Mask = ones(size(obj.InputProjections));
            else
                if FileExist(obj.filename_Mask)
                    obj.Mask = importdata(obj.filename_Mask);
                    fprintf('RESIRE: Mask used.\n');
                else
                    error('RESIRE: Mask file does not exist!')
                end
            end
            if isempty(obj.filename_InitialModel)
                fprintf('RESIRE: No initial model used.\n');
            else
                if FileExist(obj.filename_InitialModel)
                    obj.InitialModel = importdata(obj.filename_InitialModel);
                    fprintf('RESIRE: Initial model used.\n');
                else
                    error('RESIRE: Initial model file does not exist!')
                end
            end
        end
        
        function obj = CheckPrepareData(obj)
            % set number of projections
            obj.NumProjs = size(obj.InputProjections,3);
            obj.Dim1 = size(obj.InputProjections,1);
            obj.Dim2 = size(obj.InputProjections,2);
            obj.obj_dimz=round(obj.obj_dimz);
            
            % input thickness check
            if obj.obj_dimz<1
                error('RESIRE: Unrecognized sample thickness!')
            end
            
            if obj.extenedobject
                if obj.obj_dimz==1
                    [obj.obj_dimx,obj.obj_dimy,obj.obj_dimz]=objectsize(obj.Dim1,obj.Dim2,obj.Dim1,obj.vector1,obj.vector2,obj.vector3,obj.InputAngles);
                elseif obj.obj_dimz==2
                    [obj.obj_dimx,obj.obj_dimy,obj.obj_dimz]=objectsize(obj.Dim1,obj.Dim2,obj.Dim2,obj.vector1,obj.vector2,obj.vector3,obj.InputAngles);
                else
                    [obj.obj_dimx,obj.obj_dimy,obj.obj_dimz]=objectsize(obj.Dim1,obj.Dim2,obj.obj_dimz,obj.vector1,obj.vector2,obj.vector3,obj.InputAngles);
                end
            else
                obj.obj_dimx=obj.Dim1;
                obj.obj_dimy=obj.Dim2;
                if obj.obj_dimz==1
                    obj.obj_dimz=obj.Dim1;
                elseif obj.obj_dimz==2
                    obj.obj_dimz=obj.Dim2;
                end
            end
                
            % input angle and projection size check
            if obj.NumProjs~=size(obj.InputAngles,1)
                error('RESIRE: Number of projections and Angles does not match!')
            end
            
            % input angle check
            if size(obj.InputAngles,2) >3
                error('RESIRE: Input Angle 2nd dimenstion larger than three!')
            end
            
            % if only one angle set is given, make them three Euler angles
            % of y-tilt
            if size(obj.InputAngles,2) == 1
                obj.InputAngles = [zeros(obj.NumProjs ,1); obj.InputAngles; zeros(obj.NumProjs ,1)];
            end
            
            % check if interp method is legitimate
            if obj.griddingMethod ~= 1
                error('RESIRE: Unrecognized gridding method!')
            end
            
            obj.n1_oversampled = round(obj.obj_dimx * obj.oversamplingRatio_x);
            obj.n2_oversampled = round(obj.obj_dimy * obj.oversamplingRatio_y);
            obj.n3_oversampled = round(obj.obj_dimz * obj.oversamplingRatio_z);
        end
        
        function obj = runGridding(obj)
            fprintf('RESIRE: Interpolate real space projections...\n\n');
            switch obj.griddingMethod
                case 1
                    obj = obj.interp_pj_realspace();
            end
        end

        % clear big temporary arrays
        function obj = ClearCalcVariables(obj)
            obj.InputProjections = [];
            obj.InputAngles      = [];
            obj.Support          = [];
            obj.InitialModel     = [];
            obj.sum_rot_pjs      = [];
        end
        
        % save results in structure
        function SaveResults(obj)
            props = properties(obj);
            for p = 1:(numel(props)-1)
                OBJ.(props{p})=obj.(props{p});
            end
            save(obj.filename_Results, 'OBJ','-v7.3')
        end
        
        % set parameters for RESIRE class
        function obj=set_parameters(obj,varargin)
            if mod(length(varargin),2) ~= 0
                error('RESIRE: Additional argument list not divisible by 2. Options should be ''key'',''value'' pairs.')
            end
            % Apply user-provided options
            par_number = 1;
            while par_number < length(varargin)
                if isprop(obj,varargin{par_number})
                    obj.(varargin{par_number}) = varargin{par_number+1};
                else
                    error('RESIRE: Invalid option %s provided.',varargin{par_number})
                end
                par_number = par_number + 2;
            end
        end
    end
end