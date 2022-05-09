function [outlierIndexes] = findOUtNonOutliers(vector)
    %vector(isnan(vector)) = 0;
    
    nonzeroElements = ~isnan(vector);
    medianval = median(vector(nonzeroElements));
    
    % Compute the median absolute difference
    %medianval = median(vector);
    op = vector - medianval;
    outlierIndexes = or(abs(op)>0.5, isnan(op));
end