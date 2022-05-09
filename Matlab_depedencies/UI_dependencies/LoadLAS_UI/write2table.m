function retTable = write2table(lasFile)
retTable = zeros(size(lasFile.x,1),6);
retTable(:,1) = lasFile.x;
retTable(:,2) = lasFile.y;
retTable(:,3) = lasFile.z;
retTable(:,4) = get_classification(lasFile);
retTable(:,5) = lasFile.get_return_number;
retTable(:,6) = get_intensity(lasFile);
end