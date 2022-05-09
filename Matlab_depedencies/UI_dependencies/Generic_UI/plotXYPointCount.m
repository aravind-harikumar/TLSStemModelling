function plotXYPointCount(PCData,PRINT_TREE_CLOUD, GLOBAL_PRINT_ENABLE_FLAG, voxPCountArray)
    figure;
    % XZ Projection    
        subplot(3,3,1);
        plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, GLOBAL_PRINT_ENABLE_FLAG);

        totDataXZ = zeros(size(voxPCountArray,2),size(voxPCountArray,3));
        for i = 1:1:size(voxPCountArray,1);
            totDataXZ = totDataXZ +  reshape(voxPCountArray(i,:,:),[size(voxPCountArray,2) size(voxPCountArray,3)]);
        end
        subplot(3,3,2); imagesc(imrotate(totDataXZ,90)); axis square;
        totPointCnt = [];
        for ii=1:1:size(totDataXZ,1)
            totPointCnt = [totPointCnt sum(totDataXZ(ii,:))];
        end    
        totPointCnt = totPointCnt/max(totPointCnt);    
        totPointCnt(totPointCnt >0.25) = 1;
        totPointCnt(totPointCnt <=0.25) = 0;
        hold on; plot([1:1:size(totDataXZ,1)],(1-totPointCnt)*size(voxPCountArray,3),'-*','Color',[1 0 0]);

        subplot(3,3,3); 
        surf(totDataXZ,'edgecolor','none'); axis square; 
        %shading interp
    
    %YZ Projection
        subplot(3,3,4);
        plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, GLOBAL_PRINT_ENABLE_FLAG);
        %plot3(voxPCountArray(:,2),voxPCountArray(:,1),voxPCountArray(:,3),'.', 'MarkerSize' ,'8')

        totDataYZ = zeros(size(voxPCountArray,1),size(voxPCountArray,3));
        for i = 1:1:size(voxPCountArray,2);
            totDataYZ = totDataYZ +  reshape(voxPCountArray(:,i,:),[size(voxPCountArray,1) size(voxPCountArray,3)]);
        end
        subplot(3,3,5); 
        imagesc(imrotate(totDataYZ,90)); axis square;
        totPointCnt = [];
        for ii=1:1:size(totDataYZ,1)
            totPointCnt = [totPointCnt sum(totDataYZ(ii,:))];
        end
        totPointCnt = totPointCnt/max(totPointCnt);    
        totPointCnt(totPointCnt >0.25) = 1;
        totPointCnt(totPointCnt <=0.25) = 0;
        hold on; plot([1:1:size(totDataYZ,1)],(1-totPointCnt)*size(voxPCountArray,3),'-*','Color',[1 0 0]);

        subplot(3,3,6); surf(totDataYZ,'edgecolor','none'); axis square; shading interp
        %testt(totData);
    
    %XY Projection
        subplot(3,3,7);
        plotLiDARdataWithStem(PCData, PRINT_TREE_CLOUD, GLOBAL_PRINT_ENABLE_FLAG);
        totDataYZ = zeros(size(voxPCountArray,1),size(voxPCountArray,2));
        for i = 1:1:size(voxPCountArray,3);
            totDataYZ = totDataYZ +  reshape(voxPCountArray(:,:,i),[size(voxPCountArray,1) size(voxPCountArray,2)]);
        end
        subplot(3,3,8); 
        imagesc(imrotate(totDataYZ,0)); axis square;
        totPointCnt = [];
        for ii=1:1:size(totDataYZ,1)
            totPointCnt = [totPointCnt sum(totDataYZ(ii,:))];
        end
        totPointCnt = totPointCnt/max(totPointCnt);    
        totPointCnt(totPointCnt >0.25) = 1;
        totPointCnt(totPointCnt <=0.25) = 0;
        %hold on; plot([1:1:size(totDataYZ,1)],(1-totPointCnt)*size(voxPCountArray,3),'-*','Color',[1 0 0]);

        subplot(3,3,9); surf(totDataYZ,'edgecolor','none'); axis square; shading interp
        %testt(totData);
    
end

%     %totData(totData>10) = 1; totData(totData<=10) = 0;
%     axis square; hold on;
%     op1 = mat2gray(totData);
%     op = double(op1); op(op>0.5)=1;
%     imagesc(op);
%     [aa,bb]=find(op==1);
%     [centers, radii] = findCircles(op,true);
%     bestFits= EllipseDirectFit([aa bb]);
%     hold on;
%     bb = plotellipise(bestFits);
