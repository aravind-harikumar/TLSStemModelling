function StemPointArr = getStemPoints(stepSize,deg_step,PCData,vCellwise_Attrbutes,stem_center,VOXEL_SIZE, VOXEL_XYZ_DIVS,ploton)  
         StemPointArr = [];
         for i = 1:stepSize:VOXEL_XYZ_DIVS.zDiv -stepSize%;-100
            op = zeros(size(VOXEL_XYZ_DIVS.xDiv,VOXEL_XYZ_DIVS.yDiv));
            
            
            for cnt = i:1:i+stepSize
                op = op | double(mat2gray(vCellwise_Attrbutes.voxPCountArray(:,:,cnt)));
            end
            %op(op>0)=1; imagesc(op); axis equal;        
            [aa,bb]=find(op==1);
            
            %Get the closest stem point in sectors (remove outliers)
            FINAL_DATA_ARR = getsClosestPointBySector(deg_step,stem_center(i,:),[aa bb]);
            if(size(FINAL_DATA_ARR,1)>4)
                % Fit ellipse
                
                hold off;
                plot(FINAL_DATA_ARR(:,2),FINAL_DATA_ARR(:,1),'b.', 'MarkerSize', 25 ), axis equal   % Plot them
                hold on;
                plot(bb,aa,'g.', 'MarkerSize', 15 ), axis equal   % Plot them
                hold on;
                bestFits= EllipseDirectFit(FINAL_DATA_ARR); hold on;
                % plot ellipse
                [y, x] = plotellipise(bestFits,ploton);
                x = x+1; y = y-1;
                
                
                xtix = round(linspace(1,200,25)*0.02,2);
                xtixloc = linspace(1,200,25);
                set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);

                ytix = round(linspace(1,183,25)*0.02,2);
                ytixloc = linspace(1,183,25);
                set(gca,'XTickMode','auto','YTickLabel',ytix,'YTick',ytixloc);

                axis equal; axis on;
                axis([100 140 75 115]);

                xlabel('X Axis (m)','Fontname', 'Times New Roman' ,'FontSize', 20);
                ylabel('Y Axis (m)','Fontname', 'Times New Roman' ,'FontSize', 20);
                set(gca,'fontsize',20);
                
                
                % scale the image data to point cloud space
                XSpan = PCData.XMax-PCData.XMin; 
                YSpan = PCData.YMax-PCData.YMin;
                XCoord = ((x/(VOXEL_XYZ_DIVS.xDiv-1))*XSpan + VOXEL_XYZ_DIVS.xStep/2) + PCData.XMin;
                YCoord = ((y/(VOXEL_XYZ_DIVS.yDiv-1))*YSpan + VOXEL_XYZ_DIVS.yStep/2) + PCData.YMin;
                
                % store ellipse points as approximate stem points
                ZCoord = repmat(i*VOXEL_SIZE,size(x,1),1);
                StemPointArr = [StemPointArr; XCoord YCoord ZCoord];

                % prevent vertical growth og convex hull
                if(i*VOXEL_SIZE <= min(PCData.lidarDataArray(:,3)))
                   % StemPointArr(:,3) = min(PCData.lidarDataArray(:,3));
                end
                if(i*VOXEL_SIZE >= max(PCData.lidarDataArray(:,3)))
                    %StemPointArr(:,3) = max(PCData.lidarDataArray(:,3));
                end
            end
        end
    
    end