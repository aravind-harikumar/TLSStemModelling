function [stem_center,RArr,SCData,indx] = getStemCenter(PCData,HeightSlice, vCellwise_Attrbutes,VOXEL_SIZE)
    % stem_center = nan(VOXEL_XYZ_DIVS.zDiv,2);
    stem_center = []; RArr = [];
    stem_center_tmp=[]; RArrtemp=[];
    % PlotLiDARData(PCData.lidarDataArray, false, false,...
    %    PCData.htDeduction,true, 15,PCData.retMaxXYZ);
    % axis equal;
    % get scale of visualization
    xyzscale = getScaleParams(PCData,vCellwise_Attrbutes,true);

    heightstep = 100;
    for i = 1:heightstep:size(HeightSlice,3)-heightstep
        tmpImg = zeros(size(HeightSlice,1),size(HeightSlice,2));
        for j = i:1:i+heightstep  %i-(heightstep-1):1:i
            tmpImg = tmpImg + HeightSlice(:,:,j);
        end
        [centers,radii] = imfindcircles(tmpImg,[6 15],'ObjectPolarity','bright','Sensitivity',0.935);
        if(length(centers)>0)
            stem_center = [stem_center; [centers(1,1), centers(1,2), i-(heightstep/2)]];
            RArr = [RArr; radii(1)*1];
        end
    end   
    
    stem_center = stem_center(stem_center(:,3)<1300,:);
    RArr = RArr(stem_center(:,3)<1300,:);
    
    % stem center for generating cover around single-scan stems
    stem_center_tmp = stem_center;
    RArrtemp = RArr;
    CS = cat(1,0,cumsum(sqrt(sum(diff(stem_center_tmp,[],1).^2,2))));
    stem_center_tmp = interp1(CS, stem_center_tmp, unique([CS(:)' linspace(0,CS(end),length(stem_center_tmp)*30)]),'pchip');  
    RArrtemp = interp1(CS, RArrtemp, unique([CS(:)' linspace(0,CS(end),length(RArrtemp)*30)]),'pchip');  
    cylindricalvol = [];
    for i = 1:1:size(stem_center_tmp,1)
        XX = stem_center_tmp(i,1);
        YY = stem_center_tmp(i,2);
        ZZ = stem_center_tmp(i,3);
        [yunit,xunit] = circlexy1(XX,YY,RArrtemp(i)*4);
        cylindricalvol = [cylindricalvol ; [xunit' , yunit' ,repmat(ZZ,size(xunit))']];
    end    
    % get 3D point coordinates    
    stem3Dx(:,1) = stem_center_tmp(:,1)*xyzscale.ydimm + PCData.YMin;
    stem3Dx(:,2) = stem_center_tmp(:,2)*xyzscale.xdimm + PCData.XMin;
    stem3Dx(:,3) = stem_center_tmp(:,3)*xyzscale.zdimm;
    cylindricalvol(:,1) = cylindricalvol(:,1)*xyzscale.ydimm + PCData.XMin;
    cylindricalvol(:,2) = cylindricalvol(:,2)*xyzscale.xdimm + PCData.YMin;
    cylindricalvol(:,3) = cylindricalvol(:,3)*xyzscale.zdimm;
    cylindricalvol = cylindricalvol(cylindricalvol(:,3)>0,:);
    
   %if(PCData.PLOT_ON_GLOBAL_FLAG_FIG2)    
%     plot3(PCData.lidarDataArray(:,2),PCData.lidarDataArray(:,1),PCData.lidarDataArray(:,3),'*','MarkerSize', 2); axis equal; 
%     hold on;
%     plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',15,'Color',[0.5 0 0.5]);
%     hold on;
%     plot3(cylindricalvol(:,2),cylindricalvol(:,1),cylindricalvol(:,3),'*','Color',[0.5 0 0.5]); axis equal;
%     hold on;
   %plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',15,'Color',[0.5 0 0.5]);
    %end
    
    randindx = randperm(length(PCData.lidarDataArray),length(cylindricalvol));
    indx = [1+length(PCData.lidarDataArray):1:length(PCData.lidarDataArray)+length(cylindricalvol)]';
    StemCirclesArr = [PCData.lidarDataArray(randindx,1:6) [1+length(PCData.lidarDataArray):1:length(PCData.lidarDataArray)+length(cylindricalvol)]' PCData.lidarDataArray(randindx,8:13)];
    StemCirclesArr(:,1:3) = cylindricalvol;
%     PCData.lidarDataArray =  [PCData.lidarDataArray;tmp];
    SCData.lidarDataArray = StemCirclesArr;
    SCData.XMin = PCData.XMin; 
    SCData.YMin = PCData.YMin;  
    SCData.ZMin = PCData.ZMin;
    SCData.XMax = PCData.XMax;
    SCData.YMax = PCData.YMax;  
    SCData.ZMax = PCData.ZMax;
    SCData.PLOT_ON_GLOBAL_FLAG = PCData.PLOT_ON_GLOBAL_FLAG;
    
    % remove extreme points/outliers 
    stem3DxN = [];
    if(size(stem3Dx,1)>1)
        for i = 1:1:size(stem3Dx,1)
            dst = pdist2(stem3Dx(i,1:2),stem3Dx(1+1,1:2),'euclidean');
            if(dst>0.25)
               continue;
            end
            stem3DxN = [stem3DxN; stem3Dx(i,:)];
        end
    end
    %if(PCData.PLOT_ON_GLOBAL_FLAG_FIG2)
         %plot3(stem3DxN(:,1), stem3DxN(:,2), stem3DxN(:,3),'.','MarkerSize',28,'Color',[0 0.5 0]);
    %end
    
end