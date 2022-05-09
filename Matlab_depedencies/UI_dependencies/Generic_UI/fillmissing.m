function centers_fin_in = fillmissing(centers_fin_in, typeFill)
%centers_fin_in = [3, 5; 4, 5; 4 ,5; 3, 5; 3 , 6];
%centers_fin_in(3,:) = nan;
% grab integer indices of v to input to spline, and find where NaNs are
x = 1:length(centers_fin_in);
centers_fin = 0;
for colCnt = 1:1:size(centers_fin_in,2)
centers_fin = centers_fin_in(:,colCnt);
m = isnan(centers_fin);
% x(~m) contains the indices of the non-NaN values
% v(~m) contains the non-NaN values
% x(m) contains the indices of the NaN values, and thus the points at
% which we would like to query the spline interpolator
s = spline(x(~m),centers_fin(~m),x(m));
% replace NaN values with interpolated values; plot to see results
centers_fin(m) = s;
centers_fin_in(:,colCnt) = centers_fin;
end