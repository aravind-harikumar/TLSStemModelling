function plotLiDARData(lidarDataArray, gridOnOffParameterCylin,  gridOnOffParameter,htDeduction,showRawLiDARParameter,zDiv,retMaxXYZ)
if(showRawLiDARParameter==true)
    %lidarDataArray = lidarDataArray( and( lidarDataArray(:,6)>100, lidarDataArray(:,6)>150) ,:);
    plot3(lidarDataArray(:,2), lidarDataArray(:,1), lidarDataArray(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
end

if(gridOnOffParameterCylin)
    
    r=max(max(lidarDataArray(:,1:2)))+1;
    maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3));
    xamples =[];  yamples =[]; %zDiv = 15;
    for zStep = minZ :(maxZ - minZ)/zDiv : maxZ
        teta=-pi:0.314/4:pi;
        x=r*cos(teta) + retMaxXYZ(1);
        y=r*sin(teta) + retMaxXYZ(2);
        plot3(x,y,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
        plot3(x/2,y/2,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
        plot3(x/4,y/4,zeros(1,numel(x))+zStep,'Color','k', 'LineWidth',1);
        
        xamples = x(1:4:size(x,2));
        yamples = y(1:4:size(y,2));
        plot3(xamples,yamples,zeros(1,size(xamples,2))+zStep,'.','Color','k');
        
        centMat = [repmat(retMaxXYZ(:,1:2),size(xamples,2)*2,1) zeros(1,size(xamples,2)*2)'+zStep ];
        centMat([1:2:size(centMat,1)],:) = [xamples' yamples' zeros(1,size(xamples,2))'+zStep];
        
        line(centMat(:,1), centMat(:,2), centMat(:,3),'Color','k','LineWidth',1);
    end
    
    for i = 1:1:size(xamples,2)
        plot3([xamples(i) xamples(i)],[yamples(i) yamples(i)],[6.2 maxZ], '-', 'Color', 'k', 'LineWidth',1);
    end
    
end

if(gridOnOffParameter)
    
    maxX = max(lidarDataArray(:,1)); minX = min(lidarDataArray(:,1));
    maxY = max(lidarDataArray(:,2)); minY = min(lidarDataArray(:,2));
    maxZ = max(lidarDataArray(:,3)); minZ = min(lidarDataArray(:,3))-htDeduction;
    %xDiv = 10; yDiv = 10; zDiv = 25;
    
    x = minX:(maxX - minX)/zDiv: maxX;
    y = minY:(maxY - minY)/zDiv: maxY;
    z = minZ:(maxZ - minZ)/zDiv: maxZ;
    
    [X1 Y1 Z1] = meshgrid(x([1 end]),y,z);
    X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = permute(Z1,[2 1 3]);
    X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;
    
    [X2 Y2 Z2] = meshgrid(x,y([1 end]),z);
    X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;
    
    [X3 Y3 Z3] = meshgrid(x,y,z([1 end]));
    X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = permute(Z3,[3 1 2]);
    X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;
    
    %#figure('Renderer','opengl')
    h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
    set(h, 'Color',[10/256,10/256,10/256], 'LineWidth',0.5, 'LineStyle','-')
    
end
% hold off;
end

