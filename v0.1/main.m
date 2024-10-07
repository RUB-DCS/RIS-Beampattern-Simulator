clc;clear all;close all

mydir = pwd;
addpath([mydir '/optimizationFunctions'])
addpath([mydir '/RIS_pics'])
%% params

%RIS properties
param.RIS_size = 'large_bt';            % variant of RIS
param.RIS_orientation = 'standing up';  % {'laid_down'} orientation of RIS
param.N_surfaces_X = 1;                 % number of RIS modules in x-Dir
param.N_surfaces_Y = 1;                 % number of RIS modules in y-Dir   

%System properties
param.direct_channel_property = 'blocked'; % horn antenna
param.frontEdgeOffset = 3.7; % offset between horn antenna front edge and phase center in cm
param.freq = 5.53 *10^9;      % frequency, at which to optimize for

%Algorithm variables
param.opt_goal = 'maximize';
param.codebookPrecision = 1; % codebook precision in deg (~ 10/prec + 140/prec codebookPoints for localizing Tx/Rx)

param.phaseResolution = 4; % how many configurations are compared during tracking
param.UWB_angle_uncertainty = 2.5; % beam is split into +0.5*uncert deg and -0.5*uncert to cover larger area
param.UWB_dist_uncertainty = 0.1; % beam is split into +0.5*uncert m and -0.5*uncert m to cover larger area
%param.debug_initialize = 1;  % generates validation/visualization plots
param.debug_tracking = 1;    % generates validation/visualization plots

%% system knowledge to improve performance
%param.estimatedTxDist = 1.5;            % in m, estimated distance between Tx to RIS center  
% 'left'   = in front of left most RIS column or further, 
% 'center' = anywhere in front of RIS Wall
% 'right'  = in front of right most row or further, 
%param.roughTxAntennaPosition = 'left';  % w.r.t looking at RIS from front (halves the search points for Tx) 

param.Tx_dist = 2;            % in m
param.Tx_angle = -30;         % in deg
param.relative_Tx_height = 0; % zAxis offset in m (elevation w.r.t to RIS center for where Rx will be tracked)

param.relative_Rx_height = 0; % zAxis offset in m (elevation w.r.t to RIS center for where Rx will be tracked)

param.estimatedRxDist = 1.9;% distance for generating the channels for the Beamplot

trackingTechniqueInd=1; 


%choose one of those with above params
trackingTechnique{1} = 'focus';
trackingTechnique{2} = 'split';

param_focus=param;
param_split=param;



%UWB = receiver
UWB_data_position = 2;
UWB_data_angle    = 20;
%% to this in while loop of actual measurements
%for runtime = 0:10 % 10 Messungen in 10 sec

%% ToDo
% take current UWB values, this emulates movement
% UWB_data_position = UWB_data_position + rand(1,1);
% UWB_data_angle    = UWB_data_position + rand(1,1);

switch trackingTechnique {trackingTechniqueInd}
    case 'focus'
        param_focus.Rx_dist_UWB        =  UWB_data_position;         % in m
        param_focus.Rx_angle_UWB_angle =  UWB_data_angle;            % in deg
        [configSet_focus,param_focus] = trackUser(param_focus);
        param_new = param_focus;
        configSet_new = configSet_focus;

    case 'split'
        param_split.Rx_dist_UWB        =  UWB_data_position;      % in m
        param_split.Rx_angle_UWB_angle =  UWB_data_angle;         % in deg
        [configSet_split,param_split] = trackUser_BeamSplit(param_split);
        param_new = param_split;
        configSet_new = configSet_split;

end

%% visualizing calculated parameters 
%optimal RIS config
RIS_OptimalConfig = configSet_new;
figure('visible','off')
subplot(1,2,1)
imagesc((RIS_OptimalConfig>0))

%received signal at Rx
h = param_new.chan1;
g = param_new.chan2;
theta_16by16 = exp(1i*RIS_OptimalConfig);
thetaVec = reshape(theta_16by16.',[],1);

receivedSignalAtRx = pow2db( abs( h.'*diag(thetaVec)*g )^2 );
subplot(1,2,2)
bar(receivedSignalAtRx)

%generateBeampatternByHand
g_sweep = param_new.chan2_sweep;
sweepSignalStrenghts = pow2db( abs( h.'*diag(thetaVec)*g_sweep ).^2 );

Rx_sweepPoints_coor = param_new.RxSweepPoints;

%subplot(3,1,3)
fig=figure(1);
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
numSamples = 2000;
indices = round(linspace(1, length(sweepSignalStrenghts), numSamples));
sweepSignalStrenghtsSampled = sweepSignalStrenghts(indices);

[elevation_mesh, azimuth_mesh] = meshgrid(linspace(-90,90,length(sweepSignalStrenghtsSampled)), linspace(0,180,length(sweepSignalStrenghtsSampled)));

sigma_normal = 0.10;
mu_normal = 0;
values = linspace(-1,1,length(sweepSignalStrenghtsSampled));
normal = ((sigma_normal*sqrt(2*pi))^-1)*exp((-1/2)*((values-mu_normal)/sigma_normal).^2);
normal_norm = normal/max(normal);

power = sweepSignalStrenghtsSampled;
power = power + abs(min(power(:)));
power = lowpass(power, 1/(2*numSamples));
power_norm = meshgrid((power/max(power(:))).^2).*normal_norm';
power_norm = power_norm*1.7;
col_scale = power_norm.^2;

z = power_norm .* cosd(elevation_mesh) .* cosd(azimuth_mesh);
y = -power_norm .* cosd(elevation_mesh) .* sind(azimuth_mesh);
x = -power_norm .* sind(elevation_mesh);

pl5 = surf(x, y, z, col_scale);
shading flat
%shading interp; % Interpolates colors across lines and faces
colorbar;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

%scatter3(Rx_sweepPoints_coor(1,:), Rx_sweepPoints_coor(2,:), Rx_sweepPoints_coor(3,:), 50, sweepSignalStrenghts,...
%    'filled','MarkerFaceAlpha', 0.7);
%colorbar;
%end

RIS_config = (RIS_OptimalConfig>0);
%set( fig1, 'Visible', 'off' )
%% changing input of RIS
% for i=2:10
% 
% 
%     emulatedUserInput = [rand(16,16)>0.95];
%     RIS_config = xor(RIS_config ,emulatedUserInput);
%      
% 
%     thetaVec = reshape(exp(1i.*pi.*RIS_config).',[],1);
%     receivedSignalAtRx(i) = pow2db( abs( h.'*diag(thetaVec)*g )^2 );
%     sweepSignalStrenghts = pow2db( abs( h.'*diag(thetaVec)*g_sweep ).^2 );
% 
%     %visualizing
%     figure(2)
%     subplot(1,2,1)
%     imagesc(RIS_config)
%     colorbar
% 
%     subplot(1,2,2)
%     bar(receivedSignalAtRx)
% 
%     figure(1)
%     hold off
%     drawSetup(param_new)
%     scatter3(Rx_sweepPoints_coor(1,:), Rx_sweepPoints_coor(2,:), Rx_sweepPoints_coor(3,:), 50, sweepSignalStrenghts,...
%     'filled','MarkerFaceAlpha', 0.7);
%     colorbar;
% 
%     pause(0.5)
% end









