function distOfPoint2Plane =  getDistOfPoint2PlaneNew(p1,p2,p3,pext)
    retCoefVec = eqPlanenew(p1,p2,p3);
    %distOfPoint2Plane = abs(retCoefVec(1)*pext(1) + retCoefVec(2)*pext(2) + retCoefVec(3)*pext(3) + retCoefVec(4))/sqrt(sumsqr(retCoefVec(1:3)));
    dfg = size(retCoefVec,1)/4;
    dd =  reshape(retCoefVec,dfg,4);
    
    num = diag(abs(repmat([pext 1], dfg,1)*dd'));
    dem = sqrt(dd(:,1).^2 + dd(:,2).^2 + dd(:,3).^2); 
    
    distOfPoint2Plane = num./dem;
end

function retCoefVec = eqPlanenew(A,B,C)
    a = (B(:,2)-A(:,2)).*(C(:,3)-A(:,3)) - (C(:,2)-A(:,2)).*(B(:,3)-A(:,3));
    b = -((B(:,1)-A(:,1)).*(C(:,3)-A(:,3)) - (C(:,1)-A(:,1)).*(B(:,3)-A(:,3)));
    c = (B(:,1)-A(:,1)).*(C(:,2)-A(:,2)) - (C(:,1)-A(:,1)).*(B(:,2)-A(:,2));  
    d = -(a.*A(:,1) + b.*A(:,2) + c.*A(:,3));
    retCoefVec = [a;b;c;d];
end