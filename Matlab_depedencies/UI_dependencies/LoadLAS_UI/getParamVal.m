function param = getParamVal(speciesCode)
    if(strcmp(speciesCode,'ar'))
        param.POINT_DENS_THREESHOLD = 0;
        param.VOXEL_NEGH_DENSITY_THRESHOLD = 0; % - ((loopCnt-1)*4);
        param.threshDist = 0.5;
        param.MERGE_DIST = 0; 
        param.threshDistAllpoint = 0.2;
        param.MAX_ITER_RG = 10;
    elseif(strcmp(speciesCode,'la'))
        param.POINT_DENS_THREESHOLD = 0;
        param.VOXEL_NEGH_DENSITY_THRESHOLD = 0;
        param.threshDist = 0.5;
        param.MERGE_DIST = 0;
        param.threshDistAllpoint = 0.2;
        param.MAX_ITER_RG = 10;
    elseif(strcmp(speciesCode,'pc'))
        param.POINT_DENS_THREESHOLD = 0;
        param.VOXEL_NEGH_DENSITY_THRESHOLD = 0;
        param.threshDist = 0.5;
        param.MERGE_DIST = 0; 
        param.threshDistAllpoint = 0.2;
        param.MAX_ITER_RG = 10;
    end
end