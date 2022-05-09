function stemPointIndices = getTreeStemPoints(PointCluster,PCData,plotHull,plotStem)
    K = convhulln(PointCluster(:,1:3));
    if(plotStem)
        trisurf(K,PointCluster(:,2),PointCluster(:,1),PointCluster(:,3)); alpha(0.5);
    end
    % plot contained points. The output should be stem points
    stemPointIndices = inhull([PCData.lidarDataArray(:,1) PCData.lidarDataArray(:,2) PCData.lidarDataArray(:,3)] ,PointCluster(:,1:3), K);
    if(plotStem)
        plot3(PCData.lidarDataArray(stemPointIndices,2), PCData.lidarDataArray(stemPointIndices,1), PCData.lidarDataArray(stemPointIndices,3), '.','MarkerSize',5,'Color',[1 0 0]);
        %axis([-1 1 -1 1 2 4]); view(0,-90);
    end
    axis([-PCData.maxTreeWidth PCData.maxTreeWidth -PCData.maxTreeWidth PCData.maxTreeWidth 0 PCData.maxtreeHeight+(2-mod(PCData.maxtreeHeight,2));]); 
    view(-45,15);
    axis equal;
end