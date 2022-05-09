function branchPointIndx = getProminalPointsNewNew(startPoint, otherPoints, distTheshold, maxIter, branchPointIndx)
[Indx, ~] = rangesearch(otherPoints(:,1:3),startPoint(1,1:3),distTheshold);
Indx = Indx{1};
%[~,Indx] = pdist2(otherPoints(:,1:3),startPoint(1,1:3),'euclidean','Smallest',3);
testPoints = sort(otherPoints(Indx,:),5);
branchPointIndx = [branchPointIndx; testPoints(:, 4)];
if(~isempty(Indx))
    startPoint = testPoints(length(Indx),:);
    op = ismember(otherPoints(:, 4), testPoints(:, 4));
    otherPoints = removerows(otherPoints,'ind',op);

    maxIter = maxIter-1;
    try
    if((size(otherPoints(Indx,4),1) > 1) && (otherPoints(Indx(1), 5) > distTheshold) && (maxIter > 0))
        branchPointIndx = [branchPointIndx; getProminalPointsNewNew(startPoint, otherPoints, distTheshold, maxIter, branchPointIndx)];
    else
        branchPointIndx = [];
    end

    catch
        ss =0;
    end
end
end
