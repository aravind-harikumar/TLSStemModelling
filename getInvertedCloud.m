function [PCData] = getInvertedCloud(inFileFullName, OutFilePath, speciesCode,inParams)

PCData = [];
    
for file = dir(inFileFullName)'
    close all;
    fullFileName = strcat(inFileFullName,file.name); % get full file name
    if(~contains(fullFileName,"1062_TEST_TMP"))
        continue
    end
    % ----------------------------------------------------------------%
    % ----------------- LiDAR data Preprocessing -------------------- %
    % ----------------------------------------------------------------%
    % Read LiDAR data
    data = LoadData(fullFileName);
    % Preprocessing point cloud data
    PCData = dataPreProcessing_invert(data,inParams);  
    
    PCData.tmplidarDataArray = flattenProjectedCloud(PCData.tmplidarDataArray);
    % Print TLS data
%     plotLiDARdataWithStem(PCData,'asFig1','Raw TLS Data');
    VOXDiv = getVoxelDimension(PCData, 0.01);
    [csmcdmNormalized,dsmImageSegmentedArr, colorArr] =  getCSMCDM_invert(fullFileName,PCData.tmplidarDataArray, OutFilePath, VOXDiv, 2, true);
    
% 
%     surface(PCData.lidarDataArray(:,1),PCData.lidarDataArray(:,2),PCData.lidarDataArray(:,3));
%     axis equal;
end
end



%     x = min(stemcenter(:,1)):1:max(stemcenter(:,1));
%     y = min(stemcenter(:,2)):1:max(stemcenter(:,2));
%     [X,Y] = meshgrid(x,y);
%     
%     %Set the basis functions
%     f1 = ones(size(X));
%     f2 = Y.^2;
%     
%     %Write as matrix equation
%     A = [f1(:),f2(:);
%     y = stemcenter(:,3);
% %     y = data(:);
%     %Solve for coefficients
%     coeffs = A\y;
% 
% %     line = fitLine3d(stemcenter);hold on;
% % drawLine3d(line, 'color', 'm', 'LineWidth', 4);
%     CS = cat(1,0,cumsum(sqrt(sum(diff(stemcenter,[],1).^2,2))));
%     stemcenter = interp1(CS, stemcenter, unique([CS(:)' linspace(0,CS(end),100)]),'cubic');
    


