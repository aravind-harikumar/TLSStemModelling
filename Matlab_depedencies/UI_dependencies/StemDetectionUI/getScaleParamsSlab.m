function scale = getScaleParamsSlab(PCData,grid,imageto3DScale)
scale.xdimm = 0; scale.ydimm = 0;  scale.zdimm = 0;
if(imageto3DScale)
    scale.xdimm = (PCData.XMax - PCData.XMin)/(size(grid,1));
    scale.ydimm = (PCData.YMax - PCData.YMin)/(size(grid,2));
else
    scale.xdimm = size(grid,1)/(PCData.XMax - PCData.XMin);
    scale.ydimm = size(grid,2)/(PCData.YMax - PCData.YMin);
end
end