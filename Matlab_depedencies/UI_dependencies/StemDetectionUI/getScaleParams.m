function scale = getScaleParams(PCData,vCellwise_Attrbutes,imageto3DScale)
scale.xdimm = 0; scale.ydimm = 0;  scale.zdimm = 0;
if(imageto3DScale)
    scale.xdimm = (PCData.XMax - PCData.XMin)/size(vCellwise_Attrbutes.voxelCenterCellArr,1);
    scale.ydimm = (PCData.YMax - PCData.YMin)/size(vCellwise_Attrbutes.voxelCenterCellArr,2);
    scale.zdimm = (PCData.ZMax - PCData.ZMin)/size(vCellwise_Attrbutes.voxelCenterCellArr,3);
else
    scale.xdimm = size(vCellwise_Attrbutes.voxelCenterCellArr,1)/(PCData.XMax - PCData.XMin);
    scale.ydimm = size(vCellwise_Attrbutes.voxelCenterCellArr,2)/(PCData.YMax - PCData.YMin);
    scale.zdimm = size(vCellwise_Attrbutes.voxelCenterCellArr,3)/(PCData.ZMax - PCData.ZMin);
end
end