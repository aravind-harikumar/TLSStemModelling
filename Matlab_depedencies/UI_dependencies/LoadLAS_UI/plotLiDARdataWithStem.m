function plotLiDARdataWithStem(stData,inlbl,FigLabel)
if(stData.PLOT_ON_GLOBAL_FLAG)
    if(strcmp('asFig2',inlbl))
        figid = 2;
    elseif(strcmp('asFig3',inlbl))
        figid = 3;
    else
        figid = 1;
    end
    subplot(1,3,figid);
    if(stData.PRINT_STEM)
    plot3([stData.retMaxXYZ(2) stData.retMaxXYZ(2)],[stData.retMaxXYZ(1)...
        stData.retMaxXYZ(1)],[0 stData.maxtreeHeight], '-o', 'Color', [1 0 0]);
    end
    hold on;
    %To show cubic/sector grid or not ; both should not be true simulatniously
    plotLiDARData(stData.lidarDataArray, false, false,...
        stData.htDeduction,stData.PRINT_TREE_CLOUD, 15,stData.retMaxXYZ)
    
    camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); grid on;
    % maximum height
    mht = stData.maxtreeHeight + (2 - mod(stData.maxtreeHeight,2)); % round-pff to next mutiple of 2;
    % set maximum axis dimentions to be shown
    %axis([-stData.maxTreeWidth stData.maxTreeWidth -stData.maxTreeWidth ...
    % stData.maxTreeWidth 0 mht]);
    
    axis([stData.XMin-4 stData.XMax+3 stData.YMin-2 stData.YMax+2 0 ceil(stData.ZMax)]);
    
    % Set perspective and label fonts and view angle for the 3D plot
    xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
    ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
    zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);
    title(FigLabel);
    setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
    %set(gca,'XLim',[-3 3]);set(gca,'YLim',[-3 3]);set(gca,'ZLim',[0 mht]);
    %set(gca,'XTick',-3:2:3); set(gca,'YTick',-3:2:3); set(gca,'ZTick',0:2:mht)
    %set(gca,'FontSize',20)
    % set(gca,'XTickLabel',['0';'';'1';' ';'2';' ';'3';' ';'4'])
    

    
end
end