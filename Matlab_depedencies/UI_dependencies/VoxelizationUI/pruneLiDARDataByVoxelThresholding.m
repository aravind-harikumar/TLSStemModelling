function PCData = pruneLiDARDataByVoxelThresholding(PCData,voxData,VOXEL_XYZ_DIVS,threshold, pcthresh, STRATEGY, saveAsLAS,loopCnt, data)
% SelectedPointIndice = the indice of points within the selected voxels
% NearVoxlCentralPointIndice = the coordinates of the voxel centroids
% VoxelCenterPoints =  get COORDINATES of centroids of the selected voxels.

% get point indice after thesholding
%voxIndWithPoints = find( and(  voxData.voxNeighArrDensity(:)>threshold, voxData.voxPCountArray(:)> voxData.voxPCountArray(:)*pcthresh ));
voxIndWithPoints = find( and(  voxData.voxNeighArrDensity(:)>threshold, voxData.voxPCountArray(:)> pcthresh ));
SelectedPointIndice = []; NearVoxelCenterPointIndice = []; VoxelCentroids = [];
for i = 1:1:size(voxIndWithPoints)
    [iVal,jVal,kVal] = ind2sub([VOXEL_XYZ_DIVS.xDiv VOXEL_XYZ_DIVS.yDiv VOXEL_XYZ_DIVS.zDiv], voxIndWithPoints(i)); 
    SelectedPointIndice = [SelectedPointIndice; voxData.voxContainedPointIndice{iVal,jVal,kVal}];
    NearVoxelCenterPointIndice = [NearVoxelCenterPointIndice; voxData.voxNearCentroidPointIndex(iVal,jVal,kVal)]; % 
    VoxelCentroids = [VoxelCentroids; voxData.voxelCenterCellArr{iVal,jVal,kVal}];
end

if(strcmp(STRATEGY,'ALL_POINTS_IN_VOXEL'))
    indx = ismember(PCData.lidarDataArray(:,7), SelectedPointIndice);
    PCData.lidarDataArray = PCData.lidarDataArray(indx,:);  
elseif(strcmp(STRATEGY,'NEAREST_POINTS_TO_VOXEL_CENTROID'))
    indx = ismember(PCData.lidarDataArray(:,7), NearVoxelCenterPointIndice);
    PCData.lidarDataArray = PCData.lidarDataArray(indx,:);
elseif(strcmp(STRATEGY,'VOXEL_CENTROIDS_ONLY'))
    % 9 is voxel index column, and 10-12 colums are the coordinates of the voxel cetroids
    memberRowIndices = ismember(PCData.lidarDataArray(:,10:12),VoxelCentroids,'rows');
    tempLidarDataArray = unique(PCData.lidarDataArray(memberRowIndices,9:12),'rows');              
    % resizing first and then update arrays, to ovoid regenertaing unique IDS; FYI: 
    % other cols will be wrong, however useless)(temp fix)
    PCData.lidarDataArray = PCData.lidarDataArray(1:1:size(tempLidarDataArray,1),:);
    PCData.lidarDataArray(:,1:3) = tempLidarDataArray(:,2:4);
    PCData.lidarDataArray(:,9) = 0;
    PCData.lidarDataArray(:,13) = 0;
    if(loopCnt==1)
        PCData.lidarDataArrayComplete = PCData.lidarDataArray;    
    end 
else
    print('Error');
end

% build the pruned las data file.
 data = generateNewData(data,PCData.lidarDataArray);
% write .las file.
if(saveAsLAS)
     obj = write_las(data, strcat('E:\FGI_Work_LAP\2_Conifer_Species_Classification\1_Matlab_Files\1_ConiferSpeciesDetection\Conifer_classifiation_Files\LiDARDataSingleTrees\FGI_data\ar\','aa1.las'));
     fclose('all')
end

end

