
function FINAL_DATA_ARR = getsClosestPointBySector(deg_step,stem_center,X_Y)
    %stem_center = [mean(X_Y(:,1)) mean(X_Y(:,2))];
    X_Y_ANGLE = degpartition(stem_center,X_Y);    
    DIST2CEN = pdist2(stem_center,X_Y,'euclidean','Smallest',1);    
    X_Y_ANGLE_DIST2CEN = [X_Y_ANGLE DIST2CEN'];
    
    temp = removeoutliers(X_Y_ANGLE_DIST2CEN(:,3));    
    [LIA,~] = ismember(X_Y_ANGLE_DIST2CEN(:,3),temp);   
    X_Y_ANGLE_DIST2CEN = X_Y_ANGLE_DIST2CEN(LIA,:);
    
    FINAL_DATA_ARR = [];    
    for angle = 0:deg_step:360-deg_step        
       indxOfPointInSector = and(and(X_Y_ANGLE_DIST2CEN(:,3)>=angle,X_Y_ANGLE_DIST2CEN(:,3)<angle+deg_step),X_Y_ANGLE_DIST2CEN(:,4)<50);       
       X_Y_ANGLE_DIST2CEN_InSector = X_Y_ANGLE_DIST2CEN(indxOfPointInSector,:);
       [~,indsx] = min(X_Y_ANGLE_DIST2CEN_InSector(:,4));
       FINAL_DATA_ARR = [FINAL_DATA_ARR; X_Y_ANGLE_DIST2CEN_InSector(indsx,1:2)];
    end
    
%     plot(X_Y(:,1), X_Y(:,2), '*'); 
%     hold on;
%     plot(FINAL_DATA_ARR(:,1), FINAL_DATA_ARR(:,2), '*');
%     hold off;
%     
    
    
end