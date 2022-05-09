
function VDim = getVoxelDimension(stData)
    voxelSize = stData.VOXEL_SIZE;
    VDim.voxelSize = voxelSize;
    
    VDim.xDiv_const = round((stData.XMax_const - stData.XMin_const)/voxelSize,0);
    VDim.yDiv_const = round((stData.YMax_const - stData.YMin_const)/voxelSize,0);
    VDim.zDiv_const = round((stData.ZMax_const - stData.ZMin_const)/voxelSize,0);
    VDim.xStep_const = (stData.XMax_const - stData.XMin_const)/VDim.xDiv_const;
    VDim.yStep_const = (stData.XMax_const - stData.XMin_const)/VDim.yDiv_const;
    VDim.zStep_const = (stData.ZMax_const - stData.ZMin_const)/VDim.zDiv_const;
    
    VDim.xDiv = round((stData.XMax - stData.XMin)/voxelSize,0);
    VDim.yDiv = round((stData.YMax - stData.YMin)/voxelSize,0);
    VDim.zDiv = round((stData.ZMax - stData.ZMin)/voxelSize,0);
    VDim.xStep = (stData.XMax - stData.XMin)/VDim.xDiv;
    VDim.yStep = (stData.XMax - stData.XMin)/VDim.yDiv;
    VDim.zStep = (stData.ZMax - stData.ZMin)/VDim.zDiv;
    
    [VDim.rX_const,VDim.rY_const,VDim.rZ_const] = getVoxelDivisions(stData, VDim.xDiv_const, VDim.yDiv_const, VDim.zDiv_const);
    [VDim.rX,VDim.rY,VDim.rZ] = getVoxelDivisions(stData, VDim.xDiv, VDim.yDiv, VDim.zDiv);
end
