function dataAtt = dataPreProcessing_invert(singleTreeLiDARdata,inParams)

% Write normalized data to a table for performance improvement
lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
lidarDataArr(:,3) = -lidarDataArr(:,3);
minn = min(lidarDataArr(:,3));
lidarDataArr(:,3) = lidarDataArr(:,3) + minn;
%lidarDataArr = lidarDataArr(lidarDataArr(:,3)>10,:);
%lidarDataArr = [lidarDataArr;  repmat(zeros(1,size(lidarDataArr,2)), 10,1)];
%dist2center = vecnorm(lidarDataArr(:,1:2)')';
%lidarDataArr = lidarDataArr(dist2center>0.1,:);
        % for i =0:0.5:max(lidarDataArr(:,3))
        %     indx = and(lidarDataArr(:,3)>=i,lidarDataArr(:,3)<i+0.5);
        %     lidarDataArr(indx,1:2) = lidarDataArr(indx,1:2)*(1+0.2*i);   
        % end
tmplidarDataArr = lidarDataArr;
dataAtt.tmplidarDataArray = tmplidarDataArr;

if(size(lidarDataArr,1)>50000)
    lidarDataArr = lidarDataArr(randperm(size(lidarDataArr,1),50001),:);
end
treeWidth =  max(lidarDataArr(:,1)); treeWidth2 = max(lidarDataArr(:,2)); % to calculate max wisth/breadth to set plot width
dataAtt.maxTreeWidth = max(treeWidth, treeWidth2) + 0.5; % Keep it as 5 if issues arise
dataAtt.mintreeHeight = min(lidarDataArr(:,3));
dataAtt.maxtreeHeight = max(lidarDataArr(:,3));
dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;
dataAtt.htDeduction = dataAtt.treeHeight*0; % to get rid of ground noise points

% Crown height
lidarDataArr = lidarDataArr(and(lidarDataArr(:,3) > min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) < max(lidarDataArr(:,3))),:);
dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
dataAtt.maxCrownHeight = max(lidarDataArr(:,3));
dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;

%lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);
%lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);

% Identify the center point of the tree (in top view).
dataAtt.retMaxXYZ = findMaxHeightXY(lidarDataArr);

% Original LiDAR data Array
index = 1:1:size(lidarDataArr,1);
lidarDataDensityArr = ones(size(lidarDataArr,1),1); % calculateDistWeights(lidarDataArr,0.1); %
lidarDataArr = [lidarDataArr index' lidarDataDensityArr];
dataAtt.lidarDataArrayComplete =  [lidarDataArr(:,1:8) zeros(size(lidarDataArr,1),5)] ;

% Sparsified LiDAR data Array
lidarDataArrSampled = lidarDataArr(lidarDataDensityArr(:,1) >= 0.0000,:);
dataAtt.lidarDataArray =  [lidarDataArrSampled(:,1:8) zeros(size(lidarDataArrSampled,1),5 )] ;

dataAtt.XMax = max(dataAtt.lidarDataArray(:,1));
dataAtt.XMin = min(dataAtt.lidarDataArray(:,1));

dataAtt.YMax = max(dataAtt.lidarDataArray(:,2));
dataAtt.YMin = min(dataAtt.lidarDataArray(:,2));

dataAtt.ZMax = max(dataAtt.lidarDataArray(:,3));
dataAtt.ZMin = min(dataAtt.lidarDataArray(:,3));

dataAtt.XMax_const = 1;
dataAtt.XMin_const = -1;

dataAtt.YMax_const = 1;
dataAtt.YMin_const = -1;

dataAtt.ZMax_const = 10;
dataAtt.ZMin_const = 0;

dataAtt.PLOT_ON_GLOBAL_FLAG = inParams.PLOT_ON_GLOBAL_FLAG;
dataAtt.PLOT_ON_GLOBAL_FLAG_FIG2 = inParams.PLOT_ON_GLOBAL_FLAG_FIG2;
dataAtt.PRINT_STEM = inParams.PRINT_STEM; 
dataAtt.PRINT_TREE_CLOUD = inParams.PRINT_TREE_CLOUD;

dataAtt.VOXEL_SIZE = inParams.VOXEL_SIZE;

end
