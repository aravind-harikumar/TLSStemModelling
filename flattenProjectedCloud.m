function currentCloud = flattenProjectedCloud(currentCloud,In_data_name)

    currentCloud(:,3) = currentCloud(:,3) + min(currentCloud(:,3));

    
    MMS = getMaxMins(currentCloud);
    gXStep = abs(MMS.maxX - MMS.minX)/250;
    gYStep = abs(MMS.maxY - MMS.minY)/250;
    halfXGridSize = gXStep/2;       
    halfYGridSize = gXStep/2;
    %xDimGrid = size(MMS.minX+halfXGridSize : gXStep : MMS.maxX-halfXGridSize,2);
    %yDimGrid = size(MMS.minY+halfYGridSize : gYStep : MMS.maxY-halfYGridSize,2);
    %[maxVal,~] = max(currentCloud(:,3));
    minzz = MMS.minZ;
    cnt1 = 1;
    for i = MMS.minX+halfXGridSize:gXStep:MMS.maxX-halfXGridSize
        cnt2 = 1;
        for j = MMS.minY+halfYGridSize:gYStep:MMS.maxY-halfYGridSize               
            cond1 = and(currentCloud(:,1)>i-halfXGridSize,currentCloud(:,1)<=i+halfXGridSize);
            cond2 = and(currentCloud(:,2)>j-halfYGridSize,currentCloud(:,2)<=j+halfYGridSize);
            indGridPoints = find(and(cond1,cond2));
            if(size(indGridPoints,1)>0)
                localMin = min(currentCloud(indGridPoints,3));
                currentCloud(indGridPoints,3) = currentCloud(indGridPoints,3)-localMin;
            end
            cnt2 = cnt2+1;
        end
        cnt1 = cnt1+1;
    end
    
    currentCloud = currentCloud(currentCloud(:,3)<2,:);

    aa = LASread( strcat('/home/ensmingerlabgpu/Documents/MATLAB/FGI-Paper/Matlab_Code_New/LiDARDataSingleTrees/FGI_data/plot_birch/',In_data_name,'.las'));
    b = createlastemplate(aa,currentCloud);   
    LASwrite(b,strcat('/home/ensmingerlabgpu/Documents/MATLAB/FGI-Paper/Matlab_Code_New/LiDARDataSingleTrees/FGI_data/plot_birch/',In_data_name,'_b.las'));
end
