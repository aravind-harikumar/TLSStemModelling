function [csmcdmNormalized,dsmImageSegmentedArr] =  getCSMCDM(fullFileName,slabCord, OutFilePath, VOXDiv, filterOrder, plotCSM)
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
inputParams.topcutoffht = 12;
inputParams.bottomcutoffht = 3;
inputParams.slabCord = slabCord;
MMS = getMaxMins(slabCord);
inputParams.XMin = MMS.XMin;
inputParams.XMax = MMS.XMax;
inputParams.YMin = MMS.YMin;
inputParams.YMax = MMS.YMax;
inputParams.plotSurfCSM = false;
[csmcdmNormalized,dsmImageSegmentedArr] = getCSMCDMfn(fullFileName,inputParams,OutFilePath);
end

function [csm, cdm] = getCSMCDMfn(fullFileName,inP,OutFilePath)
% Initializing csm and cdm arrays
csm = zeros(inP.DimGrid1,inP.DimGrid2);
cdm = zeros(inP.DimGrid1,inP.DimGrid2);
% Generating csm and cdm
cnt1 = 1;

cond1 = and(inP.slabCord(:,3)>inP.bottomcutoffht,inP.slabCord(:,3)<inP.topcutoffht);
xyz = inP.slabCord(cond1,1:3);

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
hold on;
[centers,radii] = findCircles(raster,true);
colorbar;
axis equal;
stem_center = [];
RArr =[];
if(length(centers)>0)
    %       prev_stem_center = [centers(1,2), centers(1,1)];
    stem_center = [stem_center; [centers(:,1), centers(:,2)]];
    RArr = [RArr; radii(:)*1];
end
hold on;
plot(stem_center(:,1), stem_center(:,2),'.','MarkerSize',15,'Color',[0.5 0 0]);

xx = size(raster,2)/2;
yy = size(raster,1)/2;

subplot(1,3,1);
xyzscale = getScaleParamsSlab(inP,cdm,true);
stem3Dx(:,1) = (stem_center(:,1)-xx) * xyzscale.xdimm;
stem3Dx(:,2) = (stem_center(:,2)-yy) * xyzscale.ydimm;
plot3(stem3Dx(:,2), stem3Dx(:,1),stem3Dx(:,1)*0 + 2,'.','MarkerSize',15,'Color',[0.5 0 0]);
hold on;
axis equal;

s = LASread(fullFileName);

% crop plots
for i = 1:1:length(stem_center)
    stem3Dx(:,1) = (stem_center(:,1)-xx) * xyzscale.xdimm;
    stem3Dx(:,2) = (stem_center(:,2)-yy) * xyzscale.ydimm;
    
    % create buffer zone
    XV = stem3Dx(i,1);
    YV = stem3Dx(i,2);
    [yv,xv] = circlexy(XV,YV,1.9);
    
    % crop
    OutFile =  fullfile(OutFilePath, strcat('Spruce_2087','_',num2str(stem3Dx(i,1)),'_',num2str(stem3Dx(i,1)), '.las'));
    LASclip(s, [yv', xv'], OutFile, 'verbose', true);
end

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