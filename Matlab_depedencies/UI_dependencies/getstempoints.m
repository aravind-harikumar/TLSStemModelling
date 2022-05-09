
function [dbhSoA,DBHProposed] = getstempoints(PCDataXYZ, indx, StemCircleSlice, points, stemcenter,radius,vCellwise_Attrbutes)
%  close all;
% stemcenter = spline3dCurveInterpolation(stemcenter);
% I = points(:,:,1000);

[minval, minindx] = min(stemcenter(:,3));
stemcenter = [stemcenter; [stemcenter(minindx,1:2) -10]];

CS = cat(1,0,cumsum(sqrt(sum(diff(stemcenter,[],1).^2,2))));
stemcenter = interp1(CS, stemcenter, unique([CS(:)' linspace(0,CS(end),length(stemcenter)*1)]),'pchip');
activeindx =  find(points==1);
radius = repmat(radius(1),[size(stemcenter,1),1]);

%% Add top and bottom closing caps to stem points array.
points = addTOpBottom(points,stemcenter,radius);
ContainedPointIndice = vCellwise_Attrbutes.voxContainedPointIndice;
ContainedPointIndice = ContainedPointIndice(:);
ind = [];
for iu = 1:1:size(activeindx,1)
    ind = [ind;ContainedPointIndice{activeindx(iu)}];
end
ind = unique(ind);

posd = ~ismember(PCDataXYZ.lidarDataArray(:,7),indx);
orifdattmp = PCDataXYZ.lidarDataArray(posd,1:3);
orifdat = orifdattmp;

%% project stem to tree top
ind2 = find(points(:)>0);
[ix,iy,iz]=ind2sub(size(points),ind2);

%% Plot point cloud and stem
subplot(2,4,[1 5]);
plot3(ix,iy,iz,'.','MarkerSize',1,'Color',[0 0.5 0]); hold on;
plot3(stemcenter(:,2), stemcenter(:,1), stemcenter(:,3),'.','MarkerSize',25,'Color',[0.5 0 0]);
axis([0 size(points,2) 0 size(points,1) 0 size(points,3)]); hold on; grid on; hold on;
camproj perspective; view([-45 25]);
title('Detected Stem Center');
pbaspect([2 2 10])

%% form cylindrical seed volume around stem axis
cylindricalvol = [];
for i = 1:1:size(stemcenter,1)
    XX = stemcenter(i,1);
    YY = stemcenter(i,2);
    ZZ = stemcenter(i,3);
    [yunit,xunit] = circlexy(XX,YY,radius(i)*0.5);
%     [yunit,xunit] = circlexy(XX,YY, 3);
    cylindricalvol = [cylindricalvol ; [xunit' , yunit' ,repmat(ZZ,size(xunit))']];
end

% xmax = PCDataXYZ.XMax;
% xmin = PCDataXYZ.XMin;
% ymax = PCDataXYZ.YMax;
% ymin = PCDataXYZ.XMin;
% zmax = PCDataXYZ.ZMax;
% zmin = PCDataXYZ.ZMin;
% 
% pointa = [xmax ymax zmax;  xmin ymin zmax; xmax ymin zmax; xmin ymax zmax;
%           xmax ymax zmin;  xmin ymin zmin; xmax ymin zmin; xmin ymax zmin; ];
% xyzscale = getScaleParams(PCDataXYZ,vCellwise_Attrbutes,false);      
% pointa(:,1) = (pointa(:,1) - PCDataXYZ.XMin) * xyzscale.xdimm;
% pointa(:,2) = (pointa(:,2) - PCDataXYZ.YMin) * xyzscale.ydimm;
% pointa(:,3) = pointa(:,3)*xyzscale.zdimm;
%       
%       
% pointa = [pointa; cylindricalvol];
% points = [points; cylindricalvol];

x = cylindricalvol(:,1); y = cylindricalvol(:,2); z = cylindricalvol(:,3);
[k1,av1] = boundary(x,y,z,0.9);
%[k1,av1] = convhull(x,y,z,'simplify',true);
%trisurf(k1,x,y,z,'FaceColor','cyan'); axis([0 size(points,2) 0 size(points,1) 0 size(points,3)]); hold on; grid on; hold on;

xyzscale = getScaleParams(PCDataXYZ,vCellwise_Attrbutes,false);
stem3Dx(:,1) = (orifdat(:,1) - PCDataXYZ.XMin) * xyzscale.xdimm;
stem3Dx(:,2) = (orifdat(:,2) - PCDataXYZ.YMin) * xyzscale.ydimm;
stem3Dx(:,3) = orifdat(:,3)*xyzscale.zdimm;
% plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
% axis([0 size(points,2) 0 size(points,1) 0 size(points,3)]); hold on; grid on; hold on;

%%
% input point data
I = points;

% activeindx =  find(StemCircleSlice==1);
% I(StemCircleSlice==1) = 1;

% seed mesh data
load SphereMesh
FV.vertices = [x,y,z];
FV.faces = k1;
% subplot(1,3,3);
% title('Detected 3D Stem');
% axis([0 size(points,2) 0 size(points,1) 0 size(points,3)]); hold on; grid on; hold on;
% camproj perspective; view([-45 25]);
% pbaspect([2 2 10])

% Perform region growing
Options=struct; 
Options.Verbose=0; 
Options.Wedge=0.2;
Options.Wline=0; 
Options.Alpha=0.5; 
Options.Beta=1;
Options.Kappa=500; 
Options.Delta=1.5; 
Options.Gamma=1;
Options.Lambda = 1;
Options.Iterations=30; %Options.Sigma1=1.8; Options.Sigma2=1.8;
stemmodel1 = Snake3D(I,FV,Options);
stemmodel.vertices = stemmodel1.vertices(and(stemmodel1.vertices(:,1)>20,stemmodel1.vertices(:,2)>20),1:3);

% figure;
% tmpids = randperm(round(length(stem3Dx)/1,0));
% plot3(stem3Dx(tmpids,1), stem3Dx(tmpids,2), stem3Dx(tmpids,3),'.','MarkerSize',1,'Color',[0 0.5 0]); axis equal; grid on;
% camproj perspective; view([45 25]);  rotate3d on; axis vis3d; box on; grid on;
% xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
% ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
% zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);
% % title('TLS');
% setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
% hold on; 
% tt = stemmodel.vertices(stemmodel.vertices(:,3)<1100,:);
% plot3(tt(:,1),tt(:,2),tt(:,3),'.','MarkerSize',10); axis equal;
% % plot3(PCDataXYZ.lidarDataArray(:,1), PCDataXYZ.lidarDataArray(:,2), PCDataXYZ.lidarDataArray(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
% scalefactor = 1.3;
% stemmodelvol = ShrinkAndStreachCloud(tt, scalefactor);
% x = stemmodelvol(:,1); y = stemmodelvol(:,2); z = stemmodelvol(:,3);
% [k1,~] = boundary(x,y,z,0.7);
% trisurf(k1,x,y,z,'FaceColor','cyan'); axis([0 size(points,1) 0 size(points,2) 0 size(points,3)]); grid on;
% % plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',15,'Color',[0 0 1]);
% % camproj perspective; view([-45 25]);  alpha(0.5); axis equal;
% zoom(1.1);
% axis([-100 300 -100 300 0 2000]);
% x1 = (max(xlim)-min(xlim)); y1=(max(ylim)-min(ylim)); z1=(max(zlim)-min(zlim));
% Xlarr= 0:150:x1;  Ylarr=0:150:y1;  Zlarr=0:100:z1;
% xticks(Xlarr); yticks(Ylarr); zticks(Zlarr);
% xticklabels(num2cell( (1:1:length(Xlarr))));
% yticklabels(num2cell( (1:1:length(Ylarr))));
% zticklabels(num2cell( (1:1:length(Zlarr))));

subplot(2,4,6);
axis equal;camproj perspective; view([45 90]);
plot3(stemmodel.vertices(:,1),stemmodel.vertices(:,2),stemmodel.vertices(:,3),'*'); axis equal;

% Get DBH
subplot(2,4,7);
heightsteplow = 128; %128;
heightstepup  = 132;%132;
[dbh,dbhSoA] = getDBHatHeight(stem3Dx,stemmodel.vertices, I, heightsteplow, heightstepup, true);
dbhSoA = dbhSoA;
% axis([150 200 200 252 0 1000]);

subplot(2,4,8);
hold on;
% plot3(PCDataXYZ.lidarDataArray(:,1), PCDataXYZ.lidarDataArray(:,2), PCDataXYZ.lidarDataArray(:,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
PCDataXYZ.lidarDataArray(:,1) = (PCDataXYZ.lidarDataArray(:,1) - PCDataXYZ.XMin) * xyzscale.xdimm;
PCDataXYZ.lidarDataArray(:,2) = (PCDataXYZ.lidarDataArray(:,2) - PCDataXYZ.YMin) * xyzscale.ydimm;
PCDataXYZ.lidarDataArray(:,3) = PCDataXYZ.lidarDataArray(:,3)*xyzscale.zdimm;
pos = ~ismember(PCDataXYZ.lidarDataArray(:,7),indx);
plot3(PCDataXYZ.lidarDataArray(pos,1), PCDataXYZ.lidarDataArray(pos,2), PCDataXYZ.lidarDataArray(pos,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
axis equal;

% Scale hull and find stem points %
%[val,pos]=intersect(PCDataXYZ.lidarDataArray(:,7),indx);
% rawdata = LoadData('/home/ensmingerlabgpu/Documents/MATLAB/FGI-Paper/Matlab_Code_New/LiDARDataSingleTrees/FGI_data/Cropped_TLS/1062/Tree_1062_66.2027_-18.8505_9.las');
        subplot(1,3,1);
        scalefactor = 1;
        stemmodelvol = ShrinkAndStreachCloud(stemmodel.vertices, scalefactor);
        x = stemmodelvol(:,1); y = stemmodelvol(:,2); z = stemmodelvol(:,3);
        [k1,~] = boundary(x,y,z,0.7);
         %trisurf(k1,x,y,z,'FaceColor','cyan'); axis([0 size(points,1) 0 size(points,2) 0 size(points,3)]); hold on; grid on; hold on; axis equal;
         %plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',15,'Color',[0 0 1]);
         %alpha(0.1);
        camproj perspective; view([-45 25]);  alpha(0.5); axis equal;
        %find all points within the hull
        in1 = inhull([stem3Dx(:,1) stem3Dx(:,2) stem3Dx(:,3)],stemmodelvol);
        plot3(stem3Dx(in1,1), stem3Dx(in1,2), stem3Dx(in1,3),'.','MarkerSize',5,'Color',[0 0.5 0]); 
        axis equal;
        hold on;
%       plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',1,'Color',[0 0.5 0]);
        axis([0 300 0 300 100 800]);
        subplot(1,3,2);
        scalefactor = 1.1;
        stemmodelvol = ShrinkAndStreachCloud(stemmodel.vertices, scalefactor);
        x = stemmodelvol(:,1); y = stemmodelvol(:,2); z = stemmodelvol(:,3);
        [k1,~] = boundary(x,y,z,0.7);
         trisurf(k1,x,y,z,'FaceColor','cyan'); axis([0 size(points,1) 0 size(points,2) 0 size(points,3)]); hold on; grid on; hold on; axis equal;
        hold on;
       plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',1,'Color',[0 0.5 0]);
        camproj perspective; view([-45 25]);  alpha(0.5); axis equal;
        %find all points within the hull
        in2 = inhull([stem3Dx(:,1) stem3Dx(:,2) stem3Dx(:,3)],stemmodelvol);
        plot3(stem3Dx(in2,1), stem3Dx(in2,2), stem3Dx(in2,3),'.','MarkerSize',5,'Color',[0 0.5 0]); axis equal;
        axis([100 200 80 180 100 800]);
        sd= in1 + in2;
        subplot(1,3,3);        
       plot3(stem3Dx(sd==1,1), stem3Dx(sd==1,2), stem3Dx(sd==1,3),'.','MarkerSize',1,'Color',[0 0.5 0]);  axis equal; box on;
         view([-45 25]);
%         axis([100 200 80 180 100 800]);
         in = inhull([PCDataXYZ.lidarDataArray(:,1) PCDataXYZ.lidarDataArray(:,2) PCDataXYZ.lidarDataArray(:,3)],stemmodelvol);
         plot3(PCDataXYZ.lidarDataArray(pos,1), PCDataXYZ.lidarDataArray(pos,2), PCDataXYZ.lidarDataArray(pos,3),'.','MarkerSize',5,'Color',[0 0.5 0]);
         hold on;
        
        
        stempoints = stem3Dx(sd==1,1:3);
        stempoints = stempoints(stempoints(:,3)>180,:);
          plot3(stempoints(:,1), stempoints(:,2), stempoints(:,3),'.','MarkerSize',3,'Color',[0 0.5 0]); axis equal;
        [k1,~] = boundary(stempoints(:,1),stempoints(:,2),stempoints(:,3),0.5);
         plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',1,'Color',[0 0.5 0]); axis equal;
        hold on;
%          trisurf(k1,stempoints(:,1),stempoints(:,2),stempoints(:,3),'FaceColor','cyan'); axis([0 size(points,1) 0 size(points,2) 0 size(points,3)]); hold on; grid on; hold on; axis equal;
%          alpha(0.3);
        % oo = [stem3Dx(in,1),stem3Dx(in,2),stem3Dx(in,3)];
        hold on; grid on; hold on; axis equal;
%         plot3(stem3Dx(:,1), stem3Dx(:,2), stem3Dx(:,3),'.','MarkerSize',1,'Color',[0 0 1]); hold on;
        camproj perspective; view([45 25]);  alpha(0.5); axis equal;
%         axis([100 200 80 180 100 800]);view([45 25]); 
        scalefactor = 1;
        stempoints = ShrinkAndStreachCloud(stempoints, scalefactor);
        
        plot3(stempoints(k1,1), stempoints(k1,2), stempoints(k1,3),'.','MarkerSize',20,'Color',[0 176/256 240/256]); hold on;
%         axis([-100 400 -100 400 100 1900]);view([45 25]); 
        axis([100 200 80 180 100 800]);view([-45 25]); 
        x1 = (max(xlim)-min(xlim)); y1=(max(ylim)-min(ylim)); z1=(max(zlim)-min(zlim));
        Xlarr= 0:20:x1;  Ylarr=0:20:y1;  Zlarr=0:100:z1;
        xticks(Xlarr); yticks(Ylarr); zticks(Zlarr);
        xticklabels(num2cell( (1:1:length(Xlarr))));
        yticklabels(num2cell( (1:1:length(Ylarr))));
        zticklabels(num2cell( (1:1:length(Zlarr))));
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);



% find DBH section of the stem
poindts = stem3Dx(and(stem3Dx(:,3)>129.9,stem3Dx(:,3)<130.1),:);
% plot3(poindts(:,1),poindts(:,2),poindts(:,3),'*');
poindts = poindts(:,1:2);
camproj perspective; view([90 90]);
distarr = pdist2(poindts,poindts);
[DBHProposed,inddd] = max(distarr(:));
DBHProposed = DBHProposed*10;
title(strcat( 'DBH @ 1.3M: ', num2str(DBHProposed)));

% axis([-1 1 -1 1 1.2 1.4]);
end