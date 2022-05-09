function [csmcdmNormalized,dsmImageSegmentedArr, colorArr] =  getCSMCDM_invert(fullFileName,slabCord, OutFilePath, VOXDiv, filterOrder, plotCSM)
inputParams.DimGrid1 = size(VOXDiv.rX,2);
inputParams.inpD1 = VOXDiv.rX;
inputParams.DimGrid2 = size(VOXDiv.rY,2);
inputParams.inpD2 = VOXDiv.rY;
inputParams.gDim1Step = VOXDiv.xStep;
inputParams.gDim2Step = VOXDiv.yStep;
inputParams.halfDim1GridSize = VOXDiv.xStep*1;
inputParams.halfDim2GridSize = VOXDiv.yStep*5;
inputParams.filterOrder = filterOrder;
inputParams.plotCSM = plotCSM;
inputParams.topcutoffht = 5;
inputParams.bottomcutoffht = 4;
inputParams.slabCord = slabCord;
MMS = getMaxMins(slabCord);
inputParams.minDim1 = MMS.XMin;
inputParams.maxDim1 = MMS.XMax;
inputParams.minDim2 = MMS.YMin;
inputParams.maxDim2 = MMS.YMax;
inputParams.plotSurfCSM = false;
[csmcdmNormalized,dsmImageSegmentedArr,colorArr] = getCSMCDMfn(fullFileName,inputParams,OutFilePath);
end

function [csm, cdm, csmcdm] = getCSMCDMfn(fullFileName,inP,OutFilePath)
% Initializing csm and cdm arrays
csm = zeros(inP.DimGrid1,inP.DimGrid2);

xyz = inP.slabCord(:,1:3);



X = xyz(:,1);
Y = xyz(:,2);
[mesh, cylinders] = findTheHoles([xyz(:,1) xyz(:,2)],40);

figure();
triplot(mesh.ConnectivityList,mesh.Points(:,1),mesh.Points(:,2));

hold on;
for i = 1:length(cylinders)
    lX = X(cylinders{i});
    lY = Y(cylinders{i});
    plot([lX; lX(1)],[lY; lY(1)]);
end

dxy = 0.01;
XY = round(xyz(:,1:2) / dxy) * dxy;
xy_min = min(XY, [], 1);
xy_max = max(XY, [], 1);
xv = xy_min(1):dxy:xy_max(1);
yv = xy_min(2):dxy:xy_max(2);

[~, ~, raster] = rasterize(xyz, ...
    xv, ...
    yv, ...
    [], ...
    xyz(:,3), ...
    @(x) numel(x));

subplot(1,3,3);
imagesc(raster);


% Generating csm and cdm
cnt1 = 1;
for i = inP.minDim1+inP.halfDim1GridSize:inP.gDim1Step:inP.maxDim1-inP.halfDim1GridSize
    cnt2 = 1;
    for j = inP.minDim2+inP.halfDim2GridSize:inP.gDim2Step:inP.maxDim2-inP.halfDim2GridSize
        cond1 = and(inP.slabCord(:,1)>i-(inP.halfDim1GridSize), inP.slabCord(:,1)<=i+(inP.halfDim1GridSize));
        cond2 = and(inP.slabCord(:,2)>j-(inP.halfDim2GridSize), inP.slabCord(:,2)<=j+(inP.halfDim2GridSize));
%         if(sum(and(cond1,cond2))>0)
        csm(cnt1,cnt2) = mean(inP.slabCord(and(cond1,cond2),3));
        %imagesc(cdm); pause(0.02); hold on;
%         end
        cnt2 = cnt2+1;
    end
    cnt1 = cnt1+1;
end
csm(isnan(csm)) = 2.5;
% Intepolate CSM and CSM where values are missing
% csm = interploteImage(csm, 'linear');
% cdm = interploteImage(cdm, 'linear'); % Mean Median linear

imagesc(csm);

% Generating csm and cdm
% cnt1 = 1;

end


%    % plot the tree boundary
%    if(and(NC.ISPLOTON, inP.plotCSM))
%         f6=figure('name','CSM Plots');
%         set(f6, 'Position', [690 60 600 420]);
%         imagesc(cdm); hold on;
%         colormap(jet);
%         [heightY,lengthX] = size(csm);
%         dd = lengthX/2;
%         xTickArr = -dd-(5-mod(dd,5)):5:dd+(5-mod(dd,5));
%         yTickArr = 0:5:heightY+(5-mod(heightY,5));
%         set(gca,'XTick', xTickArr+dd+(5-mod(dd,5)));
%         set(gca,'XTickLabel', num2cell(ceil(xTickArr*inP.gDim1Step)));
%         set(gca,'YTick', yTickArr);
%         set(gca,'YTickLabel', num2cell( floor( flip(mat2gray(yTickArr)*inP.maxDim2)) ));
%         set(findall(gcf,'type','axes'),'fontsize',32);
%         set(findall(gcf,'type','text'),'fontSize',32);
%         ylabel('Tree Height','Fontname', 'Times New Roman' ,'FontSize', 36);
%         xlabel('Distance to Reference Point','Fontname', 'Times New Roman' ,'FontSize', 36);
%         title('Crown Surface Model');
%         %axis equal;
%    end