function retPointIndx = getProminalPointsNew(startPoint, otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)
[dist,Indx] = pdist2(otherPoints(:,1:3),startPoint(1,1:3),'euclidean','Smallest',NumNearestPoints);
retPointIndx = [retPointIndx; otherPoints(Indx, 5)];
otherPoints  = removerows(otherPoints,'ind',Indx(1:end-1));
maxIter = maxIter-1;
for i = 2:1:2
    if(dist(i) < distTheshold*10) && maxIter > 0
        aa =  find(otherPoints(:,5)==retPointIndx(end));
        retPointIndx = [retPointIndx;getProminalPointsNew(otherPoints(aa,:), otherPoints, distTheshold, NumNearestPoints, maxIter, retPointIndx)];
    else
        retPointIndx = [];
    end
end
end