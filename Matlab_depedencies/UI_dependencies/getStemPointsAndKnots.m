function [PCData] = getStemPointsAndKnots(PARAM)

PCData = [];

% Input TLS Path
PARAM.InputTLSPlotDataPath = strcat(PARAM.InOutBasePath,'',PARAM.InTLSDataFolder,'/');
% Cropped Tree Output path
PARAM.CroppedmergedTLSPlotPath = strcat(PARAM.InOutBasePath,'',PARAM.OutTLSDataFolder,'/',PARAM.In_data_name); % Cropped Trees Out Data Folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    CROP TLS DATA OF TREES                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(PARAM.PerformTLSCropping)
    % Read Input TLS Plot data
    PARAM.mergedTLSPlotData = LASread(strcat(PARAM.InputTLSPlotDataPath, PARAM.In_data_name,'.las'));

    % Crop TLS tree data from forest plot
    for file = dir(PARAM.inFileFullName)'
        
        close all;
        fullFileName = strcat(PARAM.inFileFullName,file.name); % get full file name
        if(~contains(fullFileName,PARAM.In_data_name))
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     TLS data Preprocessing                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read LiDAR data
        data = LoadData(fullFileName);
        % Preprocessing point cloud data
        PCData = dataPreProcessing(data,PARAM);
        % flatten point cloud
        %   PCData.tmplidarDataArray = flattenProjectedCloud(PCData.tmplidarDataArray,inParams.In_data_name);

        % generate CSM from Cloud compare
        VOXDiv = getVoxelDimension(PCData, PARAM.VOXEL_SIZE);
        %[csmcdmNormalized,dsmImageSegmentedArr] =  getCSMCDM(fullFileName,PCData.lidarDataArray, inParams.OutFilePath, VOXDiv, 2, true);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        Stem Localization                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Read the projected point cloud (done using cloudcompare for fast projection)
        ProjectedXYData = imread('/home/ensmingerlabgpu/Documents/MATLAB/FGI-Paper/Matlab_Code_New/LiDARDataSingleTrees/FGI_data/raster1029.tif');
        % normalize image
        ProjectedXYData = (ProjectedXYData-min(ProjectedXYData(:)))/(max(ProjectedXYData(:))-min(ProjectedXYData(:)));
        
        % find circular stem locations
        [centers,radii] = imfindcircles(double(ProjectedXYData),[10 44],'ObjectPolarity','bright','Sensitivity',0.91);
        distancematrix = pdist2(centers,centers,'euclidean');
        distancematrix = triu(distancematrix);
        [indxr,indxc]  = find(distancematrix>0);
        centers = centers(unique(indxr),:);
        imagesc(ProjectedXYData); hold on; plot(centers(:,1),centers(:,2),'*','Color',[0.5,0,0]);
        StemCircleAreaArray =[];
        for i =1:1:length(centers)
            [imageSizeY,imageSizeX,d] = size(ProjectedXYData);
            [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
            centerX = floor(centers(i,1));
            centerY = floor(centers(i,2));
            radius = ceil(10);
            circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
            circlePixels = mat2gray(circlePixels);
            StemCircleAreaArray = [StemCircleAreaArray; sum(sum(circlePixels.*ProjectedXYData))];
        end        
        % Remove circles with high total pixel values
        validstem = find(StemCircleAreaArray<130); % total pixel value sum inside circle > 130
        
        % imagesc(ProjectedXYData); hold on; plot(centers(ss,1),centers(ss,2),'*','Color',[0.5,0,0]);
        centers_new(:,1) = centers(validstem,1) + PCData.YMin;
        centers_new(:,2) = centers(validstem,2) + PCData.XMin;
        plot(centers(validstem,1),centers(validstem,2),'*','Color',[0.5,0,0]);
        % Normalize the data
        centers_new = [[0,0]; centers_new; [size(ProjectedXYData,1)-1,size(ProjectedXYData,2)-1]]; %add corner points for normalization
        centers_new(:,2) = size(ProjectedXYData,2) - centers_new(:,2);
        centers_new(:,1) = (centers_new(:,1)-min(centers_new(:,1)))/(max(centers_new(:,1))-min(centers_new(:,1)));
        centers_new(:,2) = (centers_new(:,2)-min(centers_new(:,2)))/(max(centers_new(:,2))-min(centers_new(:,2)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                Extract Tree-Level TLS Data                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clip the las data from the original merged TLS data
        minxx= min(PARAM.mergedTLSPlotData.record.x(:));
        maxxx= max(PARAM.mergedTLSPlotData.record.x(:));
        minyy= min(PARAM.mergedTLSPlotData.record.y(:));
        maxyy= max(PARAM.mergedTLSPlotData.record.y(:));
        centers_new(:,1) = centers_new(:,1)*(maxxx-minxx) + min(PARAM.mergedTLSPlotData.record.x(:));
        centers_new(:,2) = centers_new(:,2)*(maxyy-minyy) + min(PARAM.mergedTLSPlotData.record.y(:));
        centers_new = centers_new(2:end-1,:);
        for i =1:1:length(centers_new)
            [xv,yv] = circlexy(centers_new(i,1),centers_new(i,2),1);
            s = LASclip(PARAM.mergedTLSPlotData, [xv', yv'], fullfile(PARAM.CroppedmergedTLSPlotPath, strcat('Tree_',PARAM.In_data_name,'_', num2str(centers_new(i,1)) ,'_', num2str(centers_new(i,2)) ,'_', num2str(i),'.las')), 'verbose', true);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Perform 3D Stem Modelling and Paramter Estimation        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(PARAM.PerformStemParamterEstimation)
    outputArr = [];
    % Crop trees from plot
    for file = dir(PARAM.CroppedmergedTLSPlotPath)'        
        
        close all;
        fullFileName = fullfile(PARAM.InOutBasePath,PARAM.OutTLSDataFolder,PARAM.In_data_name, file.name); % get full file name
        % skip all non .las files
        if(~contains(fullFileName,"pine1.las"))
            continue
        end
        
        % Read Tree data
        data = LoadData(fullFileName);
        % Preprocessing point cloud data
        PCData = dataPreProcessing(data,PARAM);
        % get voxel dimensions
        VOXEL_XYZ_DIVS = getVoxelDimension(PCData,PARAM.VOXEL_SIZE);
        % get voxel attributes
        PRINT_VOXELS = true; VOXEL_TRANSPARENCY = 0.5;
        [vCellwise_Attrbutes, PCData] = formVoxelSpace(PCData,VOXEL_XYZ_DIVS,PRINT_VOXELS,VOXEL_TRANSPARENCY);
        
        % ----------------------------------------------------------------%
        %              3D stem-axis and stem-points detection             %
        % ----------------------------------------------------------------%
        HeightSliceol = vCellwise_Attrbutes.voxPCountArray;
        HeightSliceol(HeightSliceol(:)>0) = 1;
        [stem_center, RArr,StemCellData,indx] = getStemCenter(PCData,HeightSliceol,vCellwise_Attrbutes,PARAM.VOXEL_SIZE);
        
        [StemCellwise_Attrbutes, StemCellData] = formVoxelSpace(StemCellData,VOXEL_XYZ_DIVS, PRINT_VOXELS, VOXEL_TRANSPARENCY);
        StemCircleSlice = StemCellwise_Attrbutes.voxPCountArray;
        StemCircleSlice(StemCircleSlice(:)>0) = 1;
        
        vCellwise_Attrbutes1.voxelCenterCellArr = vCellwise_Attrbutes.voxelCenterCellArr;
        vCellwise_Attrbutes1.voxContainedPointIndice = vCellwise_Attrbutes.voxContainedPointIndice;
        [DBHSoA,DBHProposed] = getstempoints(PCData,indx, StemCircleSlice, HeightSliceol, stem_center, RArr, vCellwise_Attrbutes1);
        outputArr=[outputArr;  {file.name, DBHSoA, DBHProposed}];
        
    end
    
    if(exist('out1029.xls','file')==2)
        delete 'out1029.xls';
    end
    writecell(outputArr,strcat('out1029','.xls'));
    
end

end