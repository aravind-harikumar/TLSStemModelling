%calculate the relative angle for each point of the tree
function Q = degpartition(stem_center,P)    
    Q = zeros(size(P,1),3);
    Q(:,1:2) = P(:,1:2);    
    Q(:,1)= Q(:,1)-stem_center(1);
    Q(:,2)= Q(:,2)-stem_center(2);    
    dyDSM = Q(:,2); dxDSM = Q(:,1);
    % check angle (I quadrante)
    alfa=atan((dyDSM)./(dxDSM))*180/pi();
    %dist = sqrt((dyDSM.^ 2)+(dxDSM.^ 2));
    % check angle (II quadrante)
    rows = find(and(dxDSM<0,dyDSM>=0)==1);
    alfa(rows,:)=180+alfa(rows,:);
    % check angle (III quadrante)
    rows = find(and(dxDSM<0,dyDSM<0)==1);
    alfa(rows,:)=180+alfa(rows,:);
    % check angle (IV quadrante)
    rows = find(and(dxDSM>0,dyDSM<0)==1);
    alfa(rows,:)=360+alfa(rows,:);
    Q(:,3) = alfa;
    Q(:,1:2) = P(:,1:2);
    %Q(isnan(P(:,3)),3)=0;
end