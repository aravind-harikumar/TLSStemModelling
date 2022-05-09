function egffeatureArr = fitCone()
   
files = dir(strcat('.\output\','X_*.mat'));
i= 0;
for file = files'
        i = i+1;
        %subplot(6,6,i);
        hold on;
        load( fullfile('.\output\', file.name) );
        try
        X = branchPointCluster; % dataArray;    
        catch
        X = dataArray;
        end
        %Xmin = min(X(:,3));
        %X = bsxfun(@minus,X,mean(X));        
        %X(:,3) = X(:,3) + Xmin;
        plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',10,'Color',[1 0 0]);
        hold on; axis equal; grid on;
%         retEigenData = getBestEigenDirection(X);        
%         eigenVectors = retEigenData{1};
%         evalues = retEigenData{2};
%         pc1 = eigenVectors(:,1);        

%         Angle = rad2deg( atan2(norm(cross(pc1,[0 0 1])), dot(pc1,[0 0 1])) );
        %X = X*rotx(Angle);    
%         X(:,3) = X(:,3) + abs(min(X(:,3)));
%         plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',10,'Color',[1 0 0]);
%         hold on; axis equal; grid on;
%         camproj perspective; rotate3d on; axis vis3d; axis equal; axis on; view(0, 90); 
%         retEigenData = getBestEigenDirection(X);        
%         eigenVectors = retEigenData{1};
%         evalues = retEigenData{2};
%         pc1 = eigenVectors(:,1); 
%         %vectarrow([0 0 0],pc1*1);        
%         
%         width  = max( abs(min(X(:,1))-max(X(:,1))), abs(min(X(:,2))-max(X(:,2))) ) ;
%         height = max(X(:,3));
%         
        
        [n_1,~,p_1] = affine_fit(X);
        plot3(p_1(1),p_1(2),p_1(3),'ro','markersize',5,'markerfacecolor','red');
        
        [X,Y] = meshgrid(linspace(min(X(:,1)) ,max(X(:,1)),10),linspace(min(X(:,2)) ,max(X(:,2)),10));
        surf(X,Y, - (n_1(1)/n_1(3)*X+n_1(2)/n_1(3)*Y-dot(n_1,p_1)/n_1(3)),'facecolor','red','facealpha',0.5);
        
        point = p_1;
        normal = n_1';
        d = -point*normal';
        
        
        K = boundary(X(:,1),X(:,2),X(:,3),0); 
    %    trisurf(K,X(:,1),X(:,2),X(:,3));
        
        axis([-4 4 -4 4 0 30]);
%         
%         %X(:,[1 2 3]) = normalize(X(:,[1 2 3]));
%         x0 = [0; 0; height/2]; %[0 ;0; ht/2];
%         a0 = pc1; %[0;0;1];
%         phi0 = atan(width/height)*(180/pi);
%         r0 = width*2;
%         tolp = 0.001;
%         tolg = 0.001;
% 
%         X(:,3) = - X(:,3) + max( X(:,3));
%         X_Y_ofminPoint = X(find(X(:,3) == min(X(:,3)) ),1:2);
%         X(:,1) = X(:,1) - X_Y_ofminPoint(1,1)+0.01;
%         X(:,2) = X(:,2) - X_Y_ofminPoint(1,2)+0.01;
%         plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',10,'Color',[0 1 0]);
% 
% %         %X(:,1)= X(:,1)*0.3;
% %         %X(:,2)= X(:,2)*0.3;
% %         X(:,3)= -X(:,3) + max(X(:,3));
% %         X_Y_ofminPoint = X(find(X(:,3)==min(X(:,3))),1:2);
% %         X(:,1) = X(:,1) - X_Y_ofminPoint(1,1)+0.01;
% %         X(:,2) = X(:,2) - X_Y_ofminPoint(1,2)+0.01;
% %         hold on;
% %         plot3(X(:,1),X(:,2),X(:,3),'.','MarkerSize',10,'Color',[1 0 0]);
%         axis equal; grid on;
% 
%        % Cone around the z-axis, point at the origin   
%        [x0n, an, phin, rn, d, sigmah, conv, Vx0n, Van, uphin, ... 
%                 urn, GNlog, a, R0, R] = lscone(X, x0, a0, phi0, r0, ... 
%                 tolp, tolg);
% 
%         [R,A] = meshgrid(linspace(0,rn*2,110),linspace(0,2*pi,41));
%         X1 = R .* cos(A) + x0n(1);
%         Y1 = R .* sin(A) + x0n(2);
%         Z1 = R*(phin*0.5);
% 
%         X_Y_ofminPointIndx = find(Z1==min(Z1(:)));
%         X1 = X1 - X1(X_Y_ofminPointIndx(1));
%         Y1 = Y1 - Y1(X_Y_ofminPointIndx(1));
% 
%         Z1(:) = (( Z1(:) - min(Z1(:)) )/( max(Z1(:)) - min(Z1(:)) ))*height;
%         Z1 = reshape(Z1,size(X1));
% 
%         [Angle1, Angle2] = getVectorAnglesToXZandYZplanes(an);      
%         TEMPARR = [X1(:) Y1(:) Z1(:)]*rotx(Angle2);
%         X1 = reshape(TEMPARR(:,1),size(X1));
%         Y1 = reshape(TEMPARR(:,2),size(Y1));
%         Z1 = reshape(TEMPARR(:,3),size(Z1));        
%         %hold on; mesh(X1,Y1,Z1); grid on; axis equal;  
%         TEMPARR = [X1(:) Y1(:) Z1(:)]*roty(Angle1);
%         X1 = reshape(TEMPARR(:,1),size(X1));
%         Y1 = reshape(TEMPARR(:,2),size(Y1));
%         Z1 = reshape(TEMPARR(:,3),size(Z1));        
%         hold on; mesh(X1,Y1,Z1); grid on; axis equal;
%        % fet0 = (1/3)*pi*((rn)^2)*ht; % volume of cone % fet0 = (1/3)*pi*(maxx/2)^2*ht; %
   
end

end
function Norm=normalize(A)
    Norm_A=[];
    for i=1:size(A,2)
        Norm_col=(A(:,i)-min(A(:,i)))/(max(A(:,i))-min(A(:,i)));
        Norm_A=[Norm_A,Norm_col];
    end
    Norm=Norm_A;
end

function [Angle1, Angle2] = getVectorAnglesToXZandYZplanes(inputVector)  
    a = inputVector;
    nAxis = [0 0 1]; 
    xzProj = [a(1) 0 a(3)];
    yzProj = [0 a(2) a(3)];
    
    x = xzProj;
    y = nAxis; 
    Angle1 = rad2deg( 2 * atan(norm(x*norm(y) - norm(x)*y) / norm(x * norm(y) + norm(x) * y)));
    
    x = yzProj;
    y = nAxis; 
    Angle2 = rad2deg( 2 * atan(norm(x*norm(y) - norm(x)*y) / norm(x * norm(y) + norm(x) * y)));
end

function retCoefVec = eqPlanenew(A,B,C)
    a = (B(:,2)-A(:,2)).*(C(:,3)-A(:,3)) - (C(:,2)-A(:,2)).*(B(:,3)-A(:,3));
    b = -((B(:,1)-A(:,1)).*(C(:,3)-A(:,3)) - (C(:,1)-A(:,1)).*(B(:,3)-A(:,3)));
    c = (B(:,1)-A(:,1)).*(C(:,2)-A(:,2)) - (C(:,1)-A(:,1)).*(B(:,2)-A(:,2));  
    d = -(a.*A(:,1) + b.*A(:,2) + c.*A(:,3));
    retCoefVec = [a;b;c;d];
end

function distOfPoint2Plane =  getDistOfPoint2PlaneNew(p1,p2,p3,pext)
    retCoefVec = eqPlanenew(p1,p2,p3);
    %distOfPoint2Plane = abs(retCoefVec(1)*pext(1) + retCoefVec(2)*pext(2) + retCoefVec(3)*pext(3) + retCoefVec(4))/sqrt(sumsqr(retCoefVec(1:3)));
    dfg = size(retCoefVec,1)/4;
    dd =  reshape(retCoefVec,dfg,4);
    
    num = diag(abs(repmat([pext 1], dfg,1)*dd'));
    dem = sqrt(dd(:,1).^2 + dd(:,2).^2 + dd(:,3).^2); 
    
    distOfPoint2Plane = num./dem;
end

function dataAtt = dataPreProcessing(singleTreeLiDARdata)

% Write normalized data to a table for performance improvement 
        lidarDataArr = normalizeLiDARData(write2table(singleTreeLiDARdata));
        
        % Tree height
        dataAtt.maxTreeWidth = max(max(lidarDataArr(:,1)), max(lidarDataArr(:,2))); 
        dataAtt.mintreeHeight = min(lidarDataArr(:,3));
        dataAtt.maxtreeHeight = max(lidarDataArr(:,3));        
        dataAtt.treeHeight = dataAtt.maxtreeHeight - dataAtt.mintreeHeight;        
        dataAtt.htDeduction = dataAtt.treeHeight*0.1;
        
        % Crown height
        lidarDataArr = lidarDataArr(find(and(lidarDataArr(:,3) > min(lidarDataArr(:,3))+ dataAtt.htDeduction, lidarDataArr(:,3) < max(lidarDataArr(:,3)))),:);
        dataAtt.minCrownHeight = getMinCrownHeight(lidarDataArr);
        dataAtt.maxCrownHeight = max(lidarDataArr(:,3));        
        dataAtt.crownHeight = dataAtt.maxCrownHeight - dataAtt.minCrownHeight;       
        
        if(size(lidarDataArr,1)>2000)
            dataAtt.lidarDataArrSMall = lidarDataArr(randperm(size(lidarDataArr,1),2000),:);
        else
            dataAtt.lidarDataArrSMall = lidarDataArr;
        end
        
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 > 0.1),:);                
        %lidarDataArray = lidarDataArray(find(abs(lidarDataArray(:,2)).^2 + abs(lidarDataArray(:,1)).^2 < 12),:);
        
        % Identify the center point of the tree (in top view).
        dataAtt.retMaxXYZ = findMaxHeightXY(lidarDataArr);
        
        
        dataAtt.lidarDataArray = lidarDataArr;
        %dataAtt.lidarDataArray = lidarDataArray(randperm(size(lidarDataArray,1),9000),:); % for reducing the no. of samples for testing only.
        
        dataAtt.lidarDataDensityArr = zeros(size(lidarDataArr));% calculateDistWeights(lidarDataArr,30);
        dataAtt.index = 1:1:size(lidarDataArr,1);
        %dataAtt.lidarDataDensityArr = calculateDistWeights(lidarDataArr,10);

        
end

function lidarDataArray = normalizeLiDARData(lidarDataArray)
    midxPoint = (max(lidarDataArray(:,1))- min(lidarDataArray(:,1)))/2;
    midyPoint = (max(lidarDataArray(:,2))- min(lidarDataArray(:,2)))/2;
    
    lidarDataArray(:,1) = lidarDataArray(:,1)- min(lidarDataArray(:,1));
    lidarDataArray(:,2) = lidarDataArray(:,2)- min(lidarDataArray(:,2));
    lidarDataArray(:,3) = lidarDataArray(:,3)- min(lidarDataArray(:,3));

    lidarDataArray(:,1) = lidarDataArray(:,1)- midxPoint;
    lidarDataArray(:,2) = lidarDataArray(:,2)- midyPoint;
end

function crownHt = getMinCrownHeight(lidarDataArray)
   maxTreeHeight = max(lidarDataArray(:,3));
   crownHt = min(lidarDataArray(:,3));
   for i = maxTreeHeight:-1:1
       templidararray = lidarDataArray(find(and(lidarDataArray(:,3)<i,(lidarDataArray(:,3)>i-1))),1:3);
        if(size(templidararray,1)<10)
            crownHt = i;
            break;
        end
   end
end


function returnData =  LoadData(fullFileName)
    returnData = lasdata(fullFileName);
end

function retTable = write2table(lasFile)
    retTable = zeros(size(lasFile.x,1),5);
    retTable(:,1) = lasFile.x;
    retTable(:,2) = lasFile.y;
    retTable(:,3) = lasFile.z;
    retTable(:,4) = get_classification(lasFile);    
    retTable(:,5) = lasFile.get_return_number;
end

function WriteLidar2Txt(lidarDataArray,fullFileName,format)
    dlmwrite(strcat(fullFileName,format),lidarDataArray,'delimiter',' ','precision',10);
end

function isSuccess = makeLASFromTXT(exePath,exeName,inFilePath,inFileName,outFilePath,outFileName)
    command = char(strcat(exePath,exeName,{' -i '},inFilePath,inFileName,{' -o '},outFilePath,outFileName,{''}));
    [status,cmdout] = system(command); 
    isSuccess = status;
end

function retMaxXYZ = findMaxHeightXY(lidarDataArray)
    [~,indx] = max(lidarDataArray(:,3));
     retMaxXYZ = lidarDataArray(indx,1:3);
end

function retEigenData = getBestEigenDirection(inputData)
[a,~,c]=pca(inputData);
pc = zeros(3,2);

for i=1:size(a,2)
    pc(:,i) = a(:,i);%'*-c(i);
end
retEigPoints{1} = pc;
retEigPoints{2} = c;
retEigenData = retEigPoints;
end