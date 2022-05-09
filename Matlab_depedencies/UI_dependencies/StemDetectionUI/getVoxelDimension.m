function VDAtt = getVoxelDimension(inputData, voxelSize)
    VDAtt.voxelSize = voxelSize;
    VDAtt.xDiv = (inputData.XMax - inputData.XMin)/voxelSize;
    VDAtt.yDiv = (inputData.YMax - inputData.YMin)/voxelSize;
    VDAtt.zDiv = (inputData.ZMax - inputData.ZMin)/voxelSize;
    VDAtt.xStep = (inputData.XMax - inputData.XMin)/VDAtt.xDiv;
    VDAtt.yStep = (inputData.YMax - inputData.YMin)/VDAtt.yDiv;
    VDAtt.zStep = (inputData.ZMax - inputData.ZMin)/VDAtt.zDiv;
    
    % get x,y, and z voxel divisions array
    [VDAtt.rX,VDAtt.rY,VDAtt.rZ] = getVoxelDivisions(inputData, VDAtt.xDiv, VDAtt.yDiv, VDAtt.zDiv);
end