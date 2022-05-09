function [vCellData, PCData] = formVoxelSpaceStemCircles(PCData, VDim, printVoxel, voxelTransparency)

% % Plot data
% if(and(printVoxel,PCData.PLOT_ON_GLOBAL_FLAG))
%     subplot(1,3,1);
%     plotLiDARdataWithStem(PCData,'asFig2','Voxelized TLS Centroids');
% end

% %get x,y, and z voxel divisions array;
% rX = pcData.XMin + (pcData.XMax - pcData.XMin)/(xDiv-1) * [0:xDiv];
% rY = pcData.YMin + (pcData.YMax - pcData.YMin)/(yDiv-1) * [0:yDiv];
% rZ = pcData.ZMin + (pcData.ZMax - pcData.ZMin)/(zDiv-1) * [0:zDiv];

pcData.XMax 

[rX,rY,rZ] = getVoxelDivisions(PCData, VDim.xDiv-1, VDim.yDiv-1, VDim.zDiv-1);

%OT = OcTree(PCData.lidarDataArray(:,1:3));
[vCellData, PCData] = getVoxelDensityAndNebrCellIndicesNew(PCData, rX, rY, rZ, floor(VDim.xDiv), floor(VDim.yDiv), floor(VDim.zDiv));


if(and(printVoxel,PCData.PLOT_ON_GLOBAL_FLAG))
    subplot(1,3,2);
    hold on;
    %Print voxels grids;
    plot3(PCData.lidarDataArray(:,2),PCData.lidarDataArray(:,1),PCData.lidarDataArray(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
    %axis equal; camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(0, 0); %view(-45, 15);
    %axis([PCData.XMin-4 PCData.XMax+3 PCData.YMin-2 PCData.YMax+2 0 ceil(PCData.ZMax)]);

    hold on;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    cnt=1;
    for jj = 1:1:1
        %uip = uipanel('Position',[0.75 0.75 0.25 0.25]);
        %subplot(1,2,cnt)
        Indices2 = vCellData.voxNeighArrDensity(:)>=0;
        [M1,~,~] = meshgrid(rY,rX,rZ);
        [F,V,C]=ind2patch(Indices2,M1,'v'); %Creating patch data for selection of low voxels
        % Scaling and cetering
        
        YAmpFactor = (PCData.YMax - PCData.YMin)/(max(V(:,1)) - min(V(:,1)));
        XAmpFactor = (PCData.XMax - PCData.XMin)/(max(V(:,2)) - min(V(:,2)));
        ZAmpFactor = (PCData.ZMax - PCData.ZMin)/(max(V(:,3)) - min(V(:,3)));
        
        V(:,1) = V(:,1)*YAmpFactor + PCData.YMin; % + VDim.voxelSize; % origDat - shifting*scaling - centering
        V(:,2) = V(:,2)*XAmpFactor + PCData.XMin; % + VDim.voxelSize; % origDat - shifting*scaling - centering
        V(:,3) = V(:,3)*ZAmpFactor + PCData.ZMin; % + VDim.voxelSize; % origDat - shifting*scaling - centering
        
        patch('Faces',F,'Vertices',V, 'EdgeColor',[0 1 0], 'CData', ones(size(C))*100,'FaceColor',[0 0.5 0],'FaceAlpha',voxelTransparency);
        axis equal; camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; grid on; view(0, 0); %view(-45, 15);
        cnt = cnt + 1;
        axis([PCData.XMin-4 PCData.XMax+3 PCData.YMin-2 PCData.YMax+2 0 ceil(PCData.ZMax)]);
        setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
%         zoom(0.3);
    end
end