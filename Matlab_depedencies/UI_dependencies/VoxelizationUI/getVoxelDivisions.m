function[rX,rY,rZ,rX_const,rY_const,rZ_const] = getVoxelDivisions(pcData, xDiv, yDiv, zDiv)
% Perform division space using: [(minX + size_of_one_division_along_X) * scale_vector_of_size_xDiv]
rX = pcData.XMin + (pcData.XMax - pcData.XMin)/(xDiv-1) * [0:xDiv];
rY = pcData.YMin + (pcData.YMax - pcData.YMin)/(yDiv-1) * [0:yDiv];
rZ = pcData.ZMin + (pcData.ZMax - pcData.ZMin)/(zDiv-1) * [0:zDiv];

rX_const = pcData.XMin_const + (pcData.XMax_const - pcData.XMin_const)/(xDiv-1) * [0:xDiv];
rY_const = pcData.YMin_const + (pcData.YMax_const - pcData.YMin_const)/(yDiv-1) * [0:yDiv];
rZ_const = pcData.ZMin_const + (pcData.ZMax_const - pcData.ZMin_const)/(zDiv-1) * [0:zDiv];

end
