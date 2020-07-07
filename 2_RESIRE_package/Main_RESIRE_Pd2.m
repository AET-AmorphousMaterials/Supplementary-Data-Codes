%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%                     Welcome to RESIRE!                         %
%          REal Space Iterative Reconstruction Engine            %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This version incorperated extended sample recosntruction as well as 
% Radon transform for handling large samples 

addpath('./src/')
%%%%%%%%%%%%%%%%%%%%% Set RESIRE Parameters %%%%%%%%%%%%%%%%%%%%%
% See the object RESIRE_Reconstructor() for description of parameters
% please run this on super computer due to the large array size and
% oversampling ratio
RESIRE = RESIRE_reconstructor();
Path = '../1_Measured_data/Pd2_nanoparticle/';
RESIRE.filename_Projections = [Path 'projections.mat'];
RESIRE.filename_Angles = [Path 'angles.mat'];
RESIRE.filename_Support = ''; 
RESIRE.filename_InitialModel = '';
RESIRE.filename_Results = './Output/ReconObj_Pd2.mat';

RESIRE.numIterations = 200;
RESIRE.extenedobject = false;
RESIRE.obj_dimz = 1; % object thickness, 1 - same as projection dim1; 2 - same as projection dim2; Other - specified size
RESIRE.positivity = true;
RESIRE.method = 'FST'; % 'FST' or 'Radon' method for forward projection
RESIRE.oversamplingRatio_x =3.5;
RESIRE.oversamplingRatio_y =3.5;
RESIRE.oversamplingRatio_z =3.5;
RESIRE.vector1 = [0 0 1];
RESIRE.vector2 = [0 1 0];
RESIRE.vector3 = [1 0 0];
RESIRE.step_size = 2;
RESIRE.griddingMethod = 1;
RESIRE.monitor_R = 1;
RESIRE.monitorR_loopLength = 1;
RESIRE.Plot_rec = 0;
RESIRE.dtype = 'single';
%%%%%%%%%%%%%%%%%%%%% Begin RESIRE %%%%%%%%%%%%%%%%%%%%%
RESIRE = readFiles(RESIRE);
RESIRE = CheckPrepareData(RESIRE);
RESIRE = reconstruct_control(RESIRE);
Reconstruction = RESIRE.reconstruction;
save('./Output/Pd2_particle_volume.mat','Reconstruction')
SaveResults(RESIRE);