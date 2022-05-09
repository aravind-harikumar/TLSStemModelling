function fitellipsoid(xx, yy, zz)
%load('sss.mat')
% generate test data:
%xx= xx*5;
%yy= yy*5;
%zz= zz*5;
dataArray =[xx yy zz];
%[Hullk,Hullv] = convhulln(dataArray(:,1:3));

%xx = dataArray(Hullk(:,1),1);
%yy = dataArray(Hullk(:,2),2);
%zz = dataArray(Hullk(:,2),3);

%[s, t]=meshgrid([0:0.3:pi/2], [0:0.3:pi]);

% create test data:

a = (max(xx)-min(xx))/2; 
b = (max(yy)-min(yy))/2; 
c = (max(zz)-min(zz))/2; 

a=10; b=10; c=10;
% xx=a*cos(s).*cos(t);
% yy=b*cos(s).*sin(t);
% zz=c*sin(s);


% add testing noise:
%noiseIntensity = 0.3;
%xx=xx+randn(size(s))*noiseIntensity;
%yy=yy+randn(size(s))*noiseIntensity;
%zz=zz+randn(size(s))*noiseIntensity;

% do the fitting
dx=xx(:); dy=yy(:); dz=zz(:);
n=size(dx,1);

D=[dx.*dx, dy.*dy,  dz.*dz, 2.*dy.*dz, 2.*dx.*dz, 2.*dx.*dy, ...
        2.*dx, 2.*dy, 2.*dz, ones(n,1)]';

S=D*D';

v=FindFit4(S);

minX=min(dx)-1;  maxX=max(dx)+1;
minY=min(dy)-1;  maxY=max(dy)+1;
minZ=min(dz)-1;  maxZ=max(dz)+1;

% draw fitting:
nStep=200;
stepA=a/nStep; stepB=b/nStep; stepC=c/nStep;
[x, y, z]=meshgrid(minX:stepA:maxX, minY:stepB:maxY, minZ:stepC:maxZ);

try
    
    
SolidObj=v(1)*x.*x+v(2)* y.*y+v(3)*z.*z+ 2*v(4)*y.*z + 2*v(5)*x.*z + 2*v(6)*x.*y...
    + 2*v(7)*x + 2*v(8)*y + 2*v(9)*z + v(10)* ones(size(x));

%plot3(x,y,z,'.')
%plot3(dx, dy, dz, '.');

if(true)
       p = patch(isosurface(x,y,z,SolidObj, 0.0));
       isonormals(x,y,z,SolidObj, p);
       set(p, 'FaceColor', 'y', 'EdgeColor', 'none');
       daspect([1 1 1]);
       view(3);
       camlight ;
       lighting phong;

       hold on;
end
catch
    ss =0
end
end