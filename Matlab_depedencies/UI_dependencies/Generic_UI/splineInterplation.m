function splineInterplation()
    clear all; close all;
    load('matlab.mat'); load('matlab1.mat');
    
    oppAr.XMin = min(finalList(:,1));
    oppAr.XMax = max(finalList(:,1));
    oppAr.YMin = min(finalList(:,2));
    oppAr.YMax = max(finalList(:,2));
    
    VOXEL_XYZ_DIVS = getVoxelDimension(oppAr, 0.02);
    
    [cdm] = getCSMCDM(0.1, finalList*10, VOXEL_XYZ_DIVS, 2, false);
    
    x = pi*[0:.5:2];
    y = [0    0.5  0  -1   0    0.5   0;
         0.5  0   0.5  0  -0.5   0  0.5];
    pp = spline(x,y);
    yy = ppval(pp, linspace(0,2*pi,101));
    plot(yy(1,:),yy(2,:),'*b',y(1,2:5),y(2,2:5),'or')
    axis equal;
end


function [cdm] =  getCSMCDM(threshold, slabCord, VOXDiv, filterOrder, plotCSM)
    inputParams.DimGrid1 = size(VOXDiv.rX,2); inputParams.DimGrid2 = size(VOXDiv.rY,2);
    inputParams.gDim1Step = VOXDiv.xStep; inputParams.gDim2Step = VOXDiv.yStep;
    inputParams.halfDim1GridSize = VOXDiv.xStep*1; inputParams.halfDim2GridSize = VOXDiv.yStep*5;
    inputParams.filterOrder = filterOrder;
    inputParams.threshold = threshold;
    inputParams.plotCSM = plotCSM;
    inputParams.slabCord = slabCord;
    MMS = getMaxMins(slabCord);
    inputParams.minDim1 = MMS.XMin; inputParams.maxDim1 = MMS.XMax;
    inputParams.minDim2 = MMS.YMin; inputParams.maxDim2 = MMS.YMax;
    inputParams.plotSurfCSM = false;
    cdm = getCSMCDMfn(inputParams);
end


function cdm = getCSMCDMfn(inP)
    % Initializing csm and cdm arrays
    cdm = zeros(inP.DimGrid1,inP.DimGrid2);
    % Generating csm and cdm
    cnt1 = 1;
    for i = inP.minDim1+inP.halfDim1GridSize:inP.gDim1Step:inP.maxDim1-inP.halfDim1GridSize
        cnt2 = 1;
        for j = inP.minDim2+inP.halfDim2GridSize:inP.gDim2Step:inP.maxDim2-inP.halfDim2GridSize
            cond1 = and(inP.slabCord(:,1)>i-(inP.halfDim1GridSize), inP.slabCord(:,1)<=i+(inP.halfDim1GridSize));
            cond2 = and(inP.slabCord(:,2)>j-(inP.halfDim2GridSize), inP.slabCord(:,2)<=j+(inP.halfDim2GridSize));
            indGridPoints = find(and(cond1,cond2));
            if(size(indGridPoints,1)>0)
                cdm(cnt1,cnt2) = size(inP.slabCord(indGridPoints,:),1);
                %imagesc(cdm); pause(0.02); hold on;
            end
            cnt2 = cnt2+1;
        end
        cnt1 = cnt1+1;
    end

   % Intepolate CSM and CSM where values are missing
   cdm = interploteImage(cdm, 'linear'); % Mean Median linear
  
   %filterAnalysis(csm) % for test
   inP.filterOrder = 1;
   imagesc(cdm);
   cdm = imgaussfilt(cdm,inP.filterOrder);
   imagesc(cdm);
   
end

function csm = interploteImage(csm, method)
    idx = find(csm==0);
    [r,c] = ind2sub(size(csm),idx);
    for j = 1:1:size(r,1)
        avgVal = getRC(r(j),c(j),csm);
        csm(r(j),c(j)) = avgVal;
    end
    if(strcmp(method,'Mean'))
        csm(csm<mean(csm(:))) = 0; %mean(csm(:)) %0 or mean value;
    elseif(strcmp(method,'Median'))
        csm(csm<median(csm(:))) = 0; %median(csm(:)) %0 or mean value
    else
        fprintf('no thresholding performed')
    end
end

function VDAtt = getVoxelDimension(inputData, voxelSize)
    VDAtt.voxelSize = voxelSize;
    VDAtt.xDiv = round((inputData.XMax - inputData.XMin)/voxelSize,0);
    VDAtt.yDiv = round((inputData.YMax - inputData.YMin)/voxelSize,0);
    VDAtt.xStep = (inputData.XMax - inputData.XMin)/VDAtt.xDiv;
    VDAtt.yStep = (inputData.YMax - inputData.YMin)/VDAtt.yDiv;
    
    % get x,y, and z voxel divisions array
    [VDAtt.rX,VDAtt.rY] = getVoxelDivisions(inputData, VDAtt.xDiv, VDAtt.yDiv);
end

function[rX,rY] = getVoxelDivisions(pcData, xDiv, yDiv)
% Perform division space using: [(minX + size_of_one_division_along_X) * scale_vector_of_size_xDiv]
rX = pcData.XMin + (pcData.XMax - pcData.XMin)/(xDiv-1) * [0:xDiv];
rY = pcData.YMin + (pcData.YMax - pcData.YMin)/(yDiv-1) * [0:yDiv];
end


function pr = getMaxMins(slabCord)
    pr.maxX = max(slabCord(:,1));
    pr.minX = min(slabCord(:,1));
    pr.maxY = max(slabCord(:,2));
    pr.minY = min(slabCord(:,2));
    
    pr.XMax = max(slabCord(:,1));
    pr.XMin = min(slabCord(:,1));
    pr.YMax = max(slabCord(:,2));
    pr.YMin = min(slabCord(:,2));
end


function avgVal = getRC(r,c,a)
    if(and(and(r>1,r<size(a,1)),and(c>1,c<size(a,2))))
        %avgVal = (a(r,c)+a(r+1,c)+a(r-1,c)+a(r,c-1)+a(r,c+1))/4;
        avgVal = (a(r,c)+a(r+1,c)+a(r-1,c)+a(r,c-1)+a(r,c+1)+a(r+1,c+1)+a(r-1,c-1)+a(r+1,c-1)+a(r-1,c+1))/8;
    else
        avgVal = 0;
    end
end