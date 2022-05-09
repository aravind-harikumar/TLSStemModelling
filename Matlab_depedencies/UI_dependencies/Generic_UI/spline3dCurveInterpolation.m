function newCurve = spline3dCurveInterpolation(curve)
% interpote a 3d curve using spline
% path 3*x
% newPath 3*x
x = curve(:, 1);
y = curve(:, 2);
z = curve(:, 3);
t = linspace(1,10, 34)';


% Apply interpolation for each x,y and z 
tt = linspace(t(1),t(end),40);
xx = interp1(t,x,tt,'linear','extrap')';
yy = interp1(t,y,tt,'linear','extrap')';
zz = interp1(t,z,tt,'linear','extrap')';
yy = interp1(t,y,tt,'linear','extrap')';
zz = interp1(t,z,tt,'linear','extrap')';
end