function [igffeatureArr, numBranches] = getConiferIGFs_Proposed(inFileFullName, speciesCode, allFiles, printFlags)

igffeatureArr = zeros(size(allFiles,1),7);
numBranches = zeros(size(allFiles,1),1);
count = 1;

for file = allFiles'
    close all;
    fullFileName = strcat(inFileFullName,file.name); % get full file name
    
    % ----------------------------------------------------------------%
    % ----------------- LiDAR data Preprocessing -------------------- %
    % ----------------------------------------------------------------%
    % Read LiDAR data
    data = LoadData(fullFileName);
    % Preprocessing point cloud data
    PCData = dataPreProcessing(data);
    % Get parameters
    Params = getParamVal(speciesCode);
    
    % ----------------------------------------------------------------%
    %    Voxelization & cellwise attribute generation for each point  %
    % ----------------------------------------------------------------%
    subplot(1,3,1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    VOXEL_SIZE = 0.2;
    VOXEL_XYZ_DIVS = getVoxelDimension(PCData, VOXEL_SIZE);
    PRINT_VOXELS = true; VOXEL_TRANSPARENCY = 0.5;
    [vCellwise_Attrbutes, PCData] = formVoxelSpace(PCData,VOXEL_XYZ_DIVS,PRINT_VOXELS,VOXEL_TRANSPARENCY);
    PCData.lidarDataArrayComplete = PCData.lidarDataArray;    
    %Plot the conifer tree stem & also plots the complete LiDAR data
    if(printFlags.PRINTORIGINADATA)
        figure;
        subplot(1,3,2);
        PRINT_TREE_CLOUD = true;  % 1 = show LiDAR data; 0 = hide LiDAR data
        plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, printFlags.GLOBAL_PRINT_ENABLE_FLAG);
    end
    
    % ----------------------------------------------------------------%
    %                          Data pruning                           %
    % ----------------------------------------------------------------%
    % Voxel neighborhood-density thresholding based LiDAR data pruning
    PRUNE_STRATEGY = 'ALL_POINTS_IN_VOXEL'; % ['VOXEL_CENTROIDS_ONLY', 'ALL_POINTS_IN_VOXEL' and 'NEAREST_POINTS_TO_VOXEL_CENTROID']
    VOXEL_NEGH_DENSITY_THRESHOLD = 5; %Params.VOXEL_NEGH_DENSITY_THRESHOLD; 1<VOXEL_NEGH_DENSITY_THRESHOLD<8
    POINT_DENS_THREESHOLD = 0; %Params.POINT_DENS_THREESHOLD; 0<POINT_DENS_THREESHOLD<1
    SAVEASLAS = false;
    PCData = pruneLiDARDataByVoxelThresholding(PCData,vCellwise_Attrbutes, ...
        VOXEL_XYZ_DIVS, VOXEL_NEGH_DENSITY_THRESHOLD, POINT_DENS_THREESHOLD, PRUNE_STRATEGY, SAVEASLAS, 1, data);
    % Plot the conifer tree stem & also plots the complete LiDAR data
    subplot(1,3,2);
    PRINT_TREE_CLOUD = true; % 1 = show LiDAR data; 0 = hide LiDAR data
    plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, printFlags.GLOBAL_PRINT_ENABLE_FLAG);
    % link the subfigures
    ax1 = subplot(1, 3, 1); ax2 = subplot(1, 3, 2); linkSubPlotParams(ax1,ax2);
    
    % ----------------------------------------------------------------%
    %            Locates and plots the conifer branch tips            %
    % ----------------------------------------------------------------%
    plotBranchTips = true; plotConvexHull = false;
    extremeLiDARDataArray = getBoundaryPoints(PCData, speciesCode,...
        plotBranchTips, plotConvexHull, printFlags.GLOBAL_PRINT_ENABLE_FLAG);
    
    % ---------------------------------------------------------------- %
    %     This part of the code does inital region growing starting    %
    % from branche tips (but,not every point is allocated to clsuters) %
    % ---------------------------------------------------------------- %
    %load('matlab.mat');
    subplot(1,3,3);
    maxIter = Params.MAX_ITER_RG;
    distTheshold = Params.threshDist; %2.5*VOXEL_SIZE;
    [ccArray, ccIndexArray, extremeLiDARDataArray, PCData] = getConnectedComponetsNEWNEW(PCData, extremeLiDARDataArray, distTheshold, maxIter);
    ccIndexArray = sort(unique(convertIntCellDataToArray(ccIndexArray)));
    [isPointMember, MemberUniqueID, isPointNonMember, NonMemberUniqueID] = getRowIndiceWIthUniqueID(ccIndexArray, PCData.lidarDataArray,7);
    % remove points already identified to a branch using unique IDs
    PCData.lidarDataArray = PCData.lidarDataArray(isPointNonMember,:);
    
    %stemPointArr = getStemPointFromPCA(extremeLiDARDataArray, ccArray, PCData.retMaxXYZ);
    
    %load('matlab.mat');
    % get cc centroids
    subplot(1,3,3);
    %PLOT_ON = true;
    %MERGE_DIST = Params.MERGE_DIST;
    %[PCData,~] = MergeProximalcc(PCData, PLOT_ON, MERGE_DIST);
    threshDist = Params.threshDistAllpoint;
    PCData = AddRemainingPointAndPlot(PCData, threshDist);
    %ax3 = subplot(1, 3, 3); ax2 = subplot(1, 3, 2); linkSubPlotParams(ax3,ax2);
    % ----------------------------------------------------------------%
    %                          IGF Generation                         %
    % ----------------------------------------------------------------%
    figure(4);
    PRINT_TREE_CLOUD = false;  % 1 = show LiDAR data; 0 = hide LiDAR data
    plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, printFlags.GLOBAL_PRINT_ENABLE_FLAG);
    
    Hullk ={}; Hullv =[]; bc={};
    %ff = 1;
    retIGFsArray = [];
    TempArry = [PCData.lidarDataArrayComplete(:,1:3) PCData.lidarDataArrayComplete(:,13)];
    uAt = unique(TempArry(:,4));
    uAt = uAt(uAt~=0);
    
    for i=1:1:length(uAt)
        branchPointCluster = TempArry(TempArry(:,4)==uAt(i),1:3);
        upperlimit = PCData.maxtreeHeight*0.8;
        lowerlimit = PCData.maxtreeHeight*0.2;
        op = and(mean(branchPointCluster(:,3))<upperlimit, mean(branchPointCluster(:,3))>lowerlimit);
        if(and(size(branchPointCluster,1)>5, op))            
            [retFet,bcc] = getBestFitLine(branchPointCluster);
            %save( ['X_' num2str(i) '.mat'], 'branchPointCluster');            
            % retSlope, lineLength, pntDistToline, maxeval, eigRatio, totalPoints
            if(retFet.totalPoints ~= 0)
                retIGFsArray = [retIGFsArray; retFet.retSlope retFet.lineLength retFet.pntDistToline retFet.maxeval retFet.eigRatio retFet.totalPoints 0];
                hold on;
            end
        end
    end
    
    camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); grid on;
    axis([-PCData.maxTreeWidth PCData.maxTreeWidth -PCData.maxTreeWidth PCData.maxTreeWidth 0 PCData.maxtreeHeight]);
    camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(-45, 15); grid on;
    mht = PCData.maxtreeHeight + (2 - mod(PCData.maxtreeHeight,2));
    axis([-PCData.maxTreeWidth PCData.maxTreeWidth -PCData.maxTreeWidth ...
        PCData.maxTreeWidth 0 mht]);
    
    % Set perspective and label fonts and view angle for the 3D plot
    xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
    ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
    zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);
    set(gca,'XLim',[-5 5]);set(gca,'YLim',[-5 5]);set(gca,'ZLim',[0 mht]);
    set(gca,'XTick',-5:2:5); set(gca,'YTick',-5:2:5); set(gca,'ZTick',0:2:mht)
    set(gca,'FontSize',20);
    
    retIGFsVec = median(retIGFsArray,1);
    retIGFsVec(1) = retIGFsVec(1);
    retIGFsVec(2) = retIGFsVec(2)/PCData.treeHeight; % Branch Length
    retIGFsVec(4) = retIGFsVec(4)/PCData.treeHeight; % Max  eigenvalue
    retIGFsVec(6) = retIGFsVec(6)/size(PCData.lidarDataArray,1); % Total Points in a branch
    if(retIGFsVec(5)>100) % to avoid large values due to division by small number.
        retIGFsVec(5) = 100; %(eigratio)
    end
    igffeatureArr(count,:) = retIGFsVec;
    count = count + 1;
end
csvwrite( strcat(strcat(strcat('csvlist_igf_pca',num2str(count)),speciesCode ), '.csv' ),igffeatureArr);
end













































































