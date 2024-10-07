function [configSet,param] = trackUser_BeamSplit(param)
param.opt_goal = 'max_SplitBeams';
    % number of RIS Elements
    if strcmp('RIS_size','small')
        param.N_x=8;
        param.N_y=8;
    else
        param.N_x=8*2;
        param.N_y=8*2;
    end

    % position RIS in center and Rx with UWB data
   param = positionDevices(param);

    % search points for where Tx antenna might be
    %param.TxSweepPoints = generateTxSweepPoints(param);
    param.RxSplitPos = generateRxSplitPositions(param);

    % generate channels according to Tx search points
    param = generateChannels(param);

    [configsInit, gainsInit] = optimizeRIS(param);

    param.configsOpt=configsInit;

    %reshape codebook to usable format
    configSet_tracking = configsInit.quant;
%     reshapedSet = (reshape(configSet_tracking, [param.N_x*param.N_surfaces_X,param.N_y*param.N_surfaces_Y, 36]),[2,1,3]);
      if size(configSet_tracking,3)~=1;
        reshapedSet = permute(reshape(configSet_tracking, [param.N_x*param.N_surfaces_X,param.N_y*param.N_surfaces_Y, size(configSet_tracking,3)]),[2,1,3]);
      else
        reshapedSet = reshape(configSet_tracking, [param.N_x*param.N_surfaces_X,param.N_y*param.N_surfaces_Y]).';
      end
    param.codebook = reshapedSet;

    if param.debug_tracking==1
        fig2 = figure(2);
        fig2.reset;
        subplot(2,1,2)
        for k=1:size(reshapedSet,3)
            imagesc(reshapedSet(:,:,k))
        end
        pause(0.05)
    end

   configSet = reshapedSet;
    function [param] = positionDevices(param)
        % position RIS at [0,0,0]
        param.RIS_coor = [0,0,0];

        % and Rx based on UWB info

        ID_RIS_Rx =  param.Rx_dist_UWB;          % in m
        WinkRx    =  param.Rx_angle_UWB_angle;   % angle between Rx and Surface center in deg

        ID_RIS_Tx =  param.Tx_dist;          % in m
        WinkTx    =  param.Tx_angle;   % angle between Rx and Surface center in deg
        
        %include horn antenna offset
        RIS_RX_dist = ID_RIS_Rx;  
        RIS_RX_dist_offSet = RIS_RX_dist;

        RIS_TX_dist = ID_RIS_Tx + (param.frontEdgeOffset)./100;  
        RIS_TX_dist_offSet = RIS_TX_dist;

        %angular positioning relative to RIS
        Rx_X_coor = sind(WinkRx)*RIS_RX_dist_offSet;
        Rx_Y_coor = cosd(WinkRx)*RIS_RX_dist_offSet;

        Tx_X_coor = sind(WinkTx)*RIS_TX_dist_offSet;
        Tx_Y_coor = cosd(WinkTx)*RIS_TX_dist_offSet;     


        param.Rx_coor  = [Rx_X_coor,-(Rx_Y_coor),param.relative_Rx_height];
        param.Tx_coor  = [Tx_X_coor,-(Tx_Y_coor),param.relative_Tx_height];

    end
    function [TxSweepPoints] = generateTxSweepPoints(param)  
       [x,y,z] =  sampleHalfSphereEqualArea(param.estimatedTxDist);

        TxSweepPoints(1,:) = x;
        TxSweepPoints(2,:) = -y;
        TxSweepPoints(3,:)= z;
    end

    function [RxSplitPos] = generateRxSplitPositions(param)  
        
        RxSplitPos(1,1) =  sind(param.Rx_angle_UWB_angle)*(param.Rx_dist_UWB-(param.UWB_dist_uncertainty/2));
        RxSplitPos(2,1) = -cosd(param.Rx_angle_UWB_angle)*(param.Rx_dist_UWB-(param.UWB_dist_uncertainty/2)); 

        RxSplitPos(1,2) =  sind(param.Rx_angle_UWB_angle)*(param.Rx_dist_UWB+(param.UWB_dist_uncertainty/2));
        RxSplitPos(2,2) = -cosd(param.Rx_angle_UWB_angle)*(param.Rx_dist_UWB+(param.UWB_dist_uncertainty/2));

        RxSplitPos(1,3) = sind(param.Rx_angle_UWB_angle)*param.Rx_dist_UWB;
        RxSplitPos(2,3) = -cosd(param.Rx_angle_UWB_angle)*param.Rx_dist_UWB;

        RxSplitPos(3,:) = zeros(1,3);

        Pts = RxSplitPos(:,1:end-1);
        Ref = RxSplitPos(:,end);
        %dists = sqrt(sum((Pts-Ref).^2,1))
    end
    function [x, y, z] = sampleHalfSphereEqualArea(r)
        codebookPrecision = param.codebookPrecision;
        elevations = -15:codebookPrecision:-5;  % Elevation angles in degrees

        if isfield(param,'roughTxAntennaPosition')
            switch param.roughTxAntennaPosition
                case 'right'
                    azimuths = 20:codebookPrecision:20+55;  % Azimuth angles in degrees
                case 'center'
                    azimuths = 90-28:codebookPrecision:90+28;  % Azimuth angles in degreess
                case 'left'
                     azimuths = 160-55:codebookPrecision:160;  % Azimuth angles in degrees
            end
        else
            azimuths = 20:codebookPrecision:160;  % Azimuth angles in degrees
        end

        % Generate N points on a half-sphere with radius r using the "equal area" method
        [az, el] = meshgrid(deg2rad(azimuths), deg2rad(elevations));  % Generate 
        [x, y, z] = sph2cart(az, el, 1);  % Convert spherical coordinates to Cartesian coordinates
        x = x(:) * r;  % Scale points to radius r and convert x, y, z to column vectors
        y = y(:) * r;
        z = z(:) * r;
        %z = abs(z);  % Project points onto half-sphere

        %scatter3(x,-y,z)
    end

    function param = generateChannels(param)
       
        N_x = param.N_x;
        N_y = param.N_y;
        switch param.RIS_size
            case 'small'
                param.RIS_dim_x = 0.20;
                param.RIS_dim_y = 0.16;
                param.RIS_img = imread("IRS_Front.png");
            case 'large'
                param.RIS_dim_x = 0.40;
                param.RIS_dim_y = 0.32;
                param.RIS_img = imread("IRS_FrontBig.png");
            case 'large_bt'
                param.RIS_dim_x = 0.32;
                param.RIS_dim_y = 0.212;
                param.RIS_img = imread("IRS_FrontBig_new.png"); 
        end
        
        
        %% 3D virtual RIS model  

        %calculate inter-element spacings
        d_inter_x = param.RIS_dim_x/N_x;
        d_inter_y = param.RIS_dim_y/N_y;
        
        %generate inter-element positions relative to (1,1)
        xVec = repmat([1:N_x] .*d_inter_x - d_inter_x/2, N_x, 1);
        yVec = repmat([fliplr(1:N_y)] .*d_inter_y - d_inter_y/2, N_y, 1 ).';
        
        %normalize inter-element positions relative to center of RIS (0,0)
        xVec_norm = xVec - N_x/2.*d_inter_x;
        yVec_norm = yVec - N_y/2.*d_inter_y;
        
        %shift RIS center to the top left
        xShift = max(max(xVec_norm));
        yShift = max(max(yVec_norm));

        param.RealSizeX = 0.360;
        param.RealSizeY = 0.247;

        %determine space without RIS elements 
        offset_x = param.RealSizeX - param.RIS_dim_x;
        offset_y = param.RealSizeY - param.RIS_dim_y;
        
        % dynamically generate RIS wall by aligning RIS modules
        config_surfaces_x = param.N_surfaces_X;
        config_surfaces_y = param.N_surfaces_Y;
        
        temp_offset_x = 0;
        temp_offset_y = 0;
        
        xVec_norm_temp = [];
        yVec_norm_temp = [];
        
        for array_i = 1:config_surfaces_y
            for array_j = 1:config_surfaces_x
                xVec_norm_temp{array_i,array_j} = [ xVec_norm+temp_offset_x];
                yVec_norm_temp{array_i,array_j} = [ yVec_norm+temp_offset_y];
                temp_offset_x = (offset_x + param.RIS_dim_x)*array_j;
            end
            temp_offset_x = 0;
            temp_offset_y = -(offset_y + param.RIS_dim_y)*array_i;
        end

        %find new RIS center
        mid_offset_x = (config_surfaces_x-1)  * param.RealSizeX * 0.5;
        mid_offset_y = -(config_surfaces_y-1) * param.RealSizeY * 0.5;

           
        xVec_norm = cell2mat(xVec_norm_temp)-mid_offset_x;
        yVec_norm = cell2mat(yVec_norm_temp)-mid_offset_y;

        param.xVec_norm = xVec_norm;
        param.yVec_norm = yVec_norm;

        % define new position coordinates 
        switch param.RIS_orientation
            case 'laid_down'
                param.xVec_RIS = xVec_norm + param.RIS_coor(1);
                param.yVec_RIS = yVec_norm + param.RIS_coor(2);
                param.zVec_RIS = 0.*yVec_norm + param.RIS_coor(3);
            case 'standing up'
                param.xVec_RIS = xVec_norm + param.RIS_coor(1);
                param.yVec_RIS = 0.*yVec_norm + param.RIS_coor(2);
                param.zVec_RIS = yVec_norm + param.RIS_coor(3);
        end

        if param.debug_tracking == 1
            validateInitSetup(param,mid_offset_x,mid_offset_y);
        end


        Tx_coor_mat = param.Tx_coor;
        Rx_coor_mat = param.RxSplitPos.';

        xVec_RIS = reshape(param.xVec_RIS.',[],1);
        yVec_RIS = reshape(param.yVec_RIS.',[],1);
        zVec_RIS = reshape(param.zVec_RIS.',[],1);

        %calculate per element Tx-RIS and RIS-Rx distances
        xVec_RIS_Tx = xVec_RIS-Tx_coor_mat(:,1).';
        yVec_RIS_Tx = yVec_RIS-Tx_coor_mat(:,2).';
        zVec_RIS_Tx = zVec_RIS-Tx_coor_mat(:,3).';
        
        xVec_RIS_Rx = xVec_RIS-Rx_coor_mat(:,1).';
        yVec_RIS_Rx = yVec_RIS-Rx_coor_mat(:,2).';
        zVec_RIS_Rx = zVec_RIS-Rx_coor_mat(:,3).';
        
        [distanz_tr] = sqrt(xVec_RIS_Tx.^2+yVec_RIS_Tx.^2+zVec_RIS_Tx.^2);
        [distanz_rec]= sqrt(xVec_RIS_Rx.^2+yVec_RIS_Rx.^2+zVec_RIS_Rx.^2);
        
        %calculate direct distance
        distanz_dir = sqrt(sum([Rx_coor_mat-Tx_coor_mat].^2,2));
        
        
        
        %calculate channels based on distances
        c=physconst('LightSpeed');
        lambda =c/param.freq;
        
        param.chan1 = channel_near(distanz_tr,param.freq,c);  % Tx-RIS
        chan2 = channel_near(distanz_rec,param.freq,c); % RIS-Rx
        param.chan2 = chan2(:,1:end-1);
        param.chan2_ref = chan2(:,end);
        param.chanDir =  channel_near(distanz_dir,param.freq,c);
        
    end

        
 %mainAnalyticalPhaseOpt

    function fig1 = validateInitSetup(param,mid_offset_x,mid_offset_y)
    %% validate by overlaying points on RIS picture
        RIS_img = param.RIS_img;
        RIS_dim_x = param.RIS_dim_x;
        RIS_dim_y = param.RIS_dim_y;
        %get pixel dimensions of picture
        y_pix = size(RIS_img,1);
        x_pix = size(RIS_img,2);
    
        %translate picture dimensions from px to cm
        x_Ax = linspace(-RIS_dim_x/2,RIS_dim_x/2,x_pix);
        y_Ax = linspace(-RIS_dim_y/2,RIS_dim_y/2,y_pix);
    
        config_surfaces_x = param.N_surfaces_X;
        config_surfaces_y = param.N_surfaces_Y;

        offset_x = param.RealSizeX - param.RIS_dim_x;
        offset_y = param.RealSizeY - param.RIS_dim_y;
        
        % plot 3D setup 
        fig1 = figure(1);
        fig1.reset;
        grid on
        hold on

        %plot devices
        Dev_coor = [param.RIS_coor;param.Tx_coor];
        %Dev_mark{1,1} = '^';
        Dev_mark{1,1} = 'x';
        Dev_mark{2,1} = '^';
        
        for dev = 1:2
            scatter3(Dev_coor(dev,1),Dev_coor(dev,2),Dev_coor(dev,3),Dev_mark{dev,1},'LineWidth',1.5);
        end
        
        %plot RIS elements
        scatterVecX3 = param.xVec_RIS.';
        scatterVecX3 = scatterVecX3(:);
        
        scatterVecY3 = param.yVec_RIS.';
        scatterVecY3 = scatterVecY3(:);
        
        scatterVecZ3 = param.zVec_RIS.';
        scatterVecZ3 = scatterVecZ3(:);

        markerColors = ones(length(scatterVecX3), 3);
        markerColors(:,2) = linspace(1,0,length(scatterVecX3)).';
        markerColors(:,1) = linspace(1,0,length(scatterVecX3)).';
        
        pl3 = scatter3(scatterVecX3,scatterVecY3,scatterVecZ3,50,markerColors,'filled');
        pl3.LineWidth=0.01;
        pl3.MarkerEdgeColor = 'b';
        pl3.Marker = 's';
        drawnow
        scatter3(scatterVecX3(1),scatterVecY3(1),scatterVecZ3(1),'rd')
        scatter3(param.RxSplitPos(1,:),param.RxSplitPos(2,:),param.RxSplitPos(3,:),'ko')

        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Tx Localization + Beamsplit')
        legend('RIS','Rx','RIS element','1st element','Tx search point')
        
        drawnow
        view(25,25)
        axis 'equal'
    end 

   

    function [configs, gains] = optimizeRIS(params, channelParams)
        gains=[];
        switch params.opt_goal
            case 'max'
                [configs, gains] = maximizeRIS(params, channelParams);
            case 'max_CB'
                [configs] = maximizeRIS_CB(params);
            case 'max_SplitBeams'
                [configs] = maximizeRIS_SplitBeams(params);
            case 'min'
                [configs, gains] = minimizeRIS(params, channelParams);
            case 'singleElem' 
                [configs, gains] = singleElemIter(params, channelParams);
        end
        
    end
end
