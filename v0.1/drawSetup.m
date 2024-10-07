function fig1 = drawSetup(param, ax)
    %% validate by overlaying points on RIS picture

        mid_offset_x = param.mid_offset_x;
        mid_offset_y = param.mid_offset_y;

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
        if exist("ax", "var")
            fig1 = ax;
        else
            fig1 = axes();
        end
        %set(fig1,'Visible','off')
        %fig1.reset;
        grid(fig1, "on")
        hold(fig1, "on")

        %plot devices
        Dev_coor = [param.RIS_coor;param.Rx_coor];
        %Dev_mark{1,1} = '^';
        Dev_mark{1,1} = 'x';
        Dev_mark{2,1} = '^';
        
        for dev = 1:2
            scatter3(fig1, Dev_coor(dev,1),Dev_coor(dev,2),Dev_coor(dev,3),Dev_mark{dev,1},'LineWidth',1.5);
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
        
        pl3 = scatter3(fig1,scatterVecX3,scatterVecY3,scatterVecZ3,50,markerColors,'filled');
        pl3.LineWidth=0.01;
        pl3.MarkerEdgeColor = 'b';
        pl3.Marker = 's';
        drawnow
        scatter3(fig1,scatterVecX3(1),scatterVecY3(1),scatterVecZ3(1),'rd')
        
        scatter3(fig1,param.Tx_coor(1),param.Tx_coor(2),param.Tx_coor(3),'ko')

        xlabel(fig1,'x')
        ylabel(fig1,'y')
        zlabel(fig1,'z')
        title(fig1,'Beampattern')
        legend(fig1,'RIS','Rx','RIS element','1st element','Tx search point')
        
        drawnow
        view(fig1,25,25)
        axis(fig1,'equal')
    end 