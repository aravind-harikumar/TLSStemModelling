function [dbh,dbh2] = getDBHatHeight(rawdata,vertice, HeightSlice, heightsteplow, heightstepup, ploton)

vertice = vertice(vertice(:,1)~=1,:);
voxelpointsIndexInslice= find(and(vertice(:,3)>heightsteplow,vertice(:,3)<heightstepup));
[xc,yc,R,~] = circFit(vertice(voxelpointsIndexInslice,1),vertice(voxelpointsIndexInslice,2));
dbh = 2*(R*0.01)*1000;

orifdat_pointsIndexInslice = find(and(rawdata(:,3)>heightsteplow,rawdata(:,3)<heightstepup));
rawsubset = [rawdata(orifdat_pointsIndexInslice,1), rawdata(orifdat_pointsIndexInslice,2)];
ampFact  =1.3;
vertice1 = vertice;
[k1,av1] = boundary(vertice1(:,1),vertice1(:,2),0.5);
x= vertice1(:,1); 
y = vertice1(:,2);
plot(x(k1),y(k1),'*r');
hold on;
axis([10 400 0 400]);

vertice2 = vertice1*ampFact;
[k11,av1] = boundary(vertice2(:,1),vertice2(:,2),0.5);
diffx = abs( min(vertice2(:,1)) - min(vertice1(:,1)) );
diffy = abs(min(vertice2(:,2)) - min(vertice1(:,2)));
vertice2(:,1) = vertice2(:,1) - diffx; 
vertice2(:,2) = vertice2(:,2) - diffy;

mm1 = (max(vertice2(:,1)) - max(vertice1(:,1)))/2;
mm2 = (max(vertice2(:,2)) - max(vertice1(:,2)))/2;

vertice2(:,1)  = vertice2(:,1) -mm1;
vertice2(:,2) = vertice2(:,2)  -mm2;

x2= vertice2(:,1); 
y2 = vertice2(:,2);
plot(x2(k11),y2(k11),'*b');
axis([10 400 0 400]);
axis equal;
hold on; grid on; hold on;

newBoun = [x2(k11),y2(k11)];
in = inpolygon(rawsubset(:,1),rawsubset(:,2),newBoun(:,1),newBoun(:,2));
hold on;
plot(rawsubset(in,1),rawsubset(in,2),'.g')


[xc1,yc1,R1,~] = circFit(rawsubset(in,1),rawsubset(in,2));
dbh1 = 2*(R1*0.01)*1000;


[xc1,yc1,R2,~] = circFit(rawsubset(in,1),rawsubset(in,2));
dbh2 = 2*(R2*0.01)*1000;

if(ploton)
%     figure;
    tmpImg = zeros(size(HeightSlice,1),size(HeightSlice,2));
    for j = heightsteplow:1:heightstepup
        tmpImg = tmpImg + HeightSlice(:,:,j);
    end
        %     subplot(2,4,6);
        %     imagesc(tmpImg)
        %     c = [yc,xc];
        %     hold on;
        %     viscircles(c,R);
        %     hold on;
        %     plot3(vertice(voxelpointsIndexInslice,2),vertice(voxelpointsIndexInslice,1),vertice(voxelpointsIndexInslice,3),'.')
        %     axis equal;
        %     title(strcat('Circle Fit on Convex Hull Vertices @ 1.3M: ', num2str(dbh)))
        %     xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        %     ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        % %     zlabel('Tree height','Fontname', 'Times New Roman' ,'FontSize', 14);
        % %     setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
        %     axis([0 300 0 300 0 1000]);
    
        subplot(2,4,7);
        c1 = [yc1,xc1];
        imagesc(zeros(size(tmpImg)))
        hold on;
        viscircles(c1,R1);
        hold on;
        plot(rawsubset(in,2),rawsubset(in,1),'.g')
        axis equal;
        title(strcat( 'Circle Fit on Raw Data Points @ 1.3M: ', num2str(dbh1)))
        xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
        axis([0 size(tmpImg,1) 0 size(tmpImg,2)]);
    
%     subplot(2,4,7);
%     c1 = [yc1,xc1];
%     imagesc(zeros(size(tmpImg)))
%     hold on;
%     viscircles(c1,R1);
%     hold on;
%     plot(rawsubset(in,2),rawsubset(in,1),'.g')
%     axis equal;
%     title(strcat( 'Circle Fit on Raw Data Points @ 1.3M: ', num2str(dbh2)))
%     xlabel('X Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
%     ylabel('Y Axis','Fontname', 'Times New Roman' ,'FontSize', 14);
%     axis([0 300 0 300]);
end
end

