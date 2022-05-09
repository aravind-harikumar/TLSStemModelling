
function newpointarr = ShrinkAndStreachCloud(PointCluster,scalefactor)
%      load('matlabttt.mat')
%      scalefactor = 2;
%      
    tmpdd = PointCluster;
%     plot3(tmpdd(:,1), tmpdd(:,2), tmpdd(:,3), '.','MarkerSize',5);  axis equal;  
    
    newpointarr = [];
    for ii=0:50:max(tmpdd(:,3))
    
    PointCluster = tmpdd(and(tmpdd(:,3)>ii,tmpdd(:,3)<=ii+100),:);    
    if(length(PointCluster)<4)
        continue;
    end
    
    tmp = PointCluster;
    PointCluster1= [];
    PointCluster1(:,1:2) = PointCluster(:,1:2)*scalefactor;
    PointCluster1(:,3) = PointCluster(:,3);

    [maxpoint,~] = max(tmp); % MAX Point A
    [scaledmaxpoint,~] = max(PointCluster1); % MAX Point A (in scaled)

    [minpoint,~] = min(tmp); % MIN Point B
    [scaledminpoint,~] = min(PointCluster1); % MIN Point B (in scaled)%     
%     hold on;
     
    shiftx = scaledminpoint(1) - minpoint(1);
    shifty = scaledminpoint(2) - minpoint(2); 
    
    shidtedScaledCloud=[];
    shidtedScaledCloud(:,1) = PointCluster1(:,1) - shiftx;
    shidtedScaledCloud(:,2) = PointCluster1(:,2) - shifty;
      
    [maxpointSS,~] = max(shidtedScaledCloud(:,1:2));
    [minpointSS,~] = min(shidtedScaledCloud(:,1:2)); 
  
    xratio = (maxpointSS(1)-maxpoint(1));
    yratio = (maxpointSS(2)-maxpoint(2));
    
%     hold on;
%     plot(maxpointSS(:,1), maxpointSS(:,2), '.','MarkerSize',15,'Color',[1 0 0]);  axis equal; 
%     plot(minpointSS(:,1), minpointSS(:,2), '.','MarkerSize',15,'Color',[0 1 0]);  axis equal; 
%     
%     plot(maxpoint(:,1), maxpoint(:,2), '.','MarkerSize',15,'Color',[0 0 1]);  axis equal; 
%     plot(minpoint(:,1), minpoint(:,2), '.','MarkerSize',45,'Color',[0 0.6 0.2]);  axis equal; 
%     
    PointCluster2 =[];
    PointCluster2(:,1) = shidtedScaledCloud(:,1) - xratio/2;
    PointCluster2(:,2) = shidtedScaledCloud(:,2) - yratio/2;
    PointCluster2(:,3) = PointCluster1(:,3); 
    
    newpointarr = [newpointarr;PointCluster2];
    
    end
    
%     hold on;
%     plot3(newpointarr(:,1), newpointarr(:,2), newpointarr(:,3), '.','MarkerSize',5,'Color',[1 0 0]);  axis equal;  
%     
%     xdiff = (max(PointCluster1(:,1)) - max(tmp(:,1))) - (max(PointCluster1(:,1))-min(PointCluster1(:,1)))*0.5 + (max(tmp(:,1))-min(tmp(:,1)))*0.5;
%     ydiff = (max(PointCluster1(:,2)) - max(tmp(:,2))) - (max(PointCluster1(:,2))-min(PointCluster1(:,2)))*0.5 + (max(tmp(:,2))-min(tmp(:,2)))*0.5;
%     %zdiff = (max(PointCluster(:,3)) - max(tmp(:,3))) - (max(PointCluster(:,3))-min(PointCluster(:,3)))*0.5 + (max(tmp(:,3))-min(tmp(:,3)))*0.5;      
%     PointCluster2(:,1) = PointCluster1(:,1) - ydiff;
%     PointCluster2(:,2) = PointCluster1(:,2) - xdiff;
%     PointCluster2(:,3) = PointCluster1(:,3);
%     hold on;
%      plot3(PointCluster2(:,1), PointCluster2(:,2), PointCluster2(:,3), '.','MarkerSize',5,'Color',[1 0 0]);  axis equal;  
end