function points = addTOpBottom(points,stemcenter,RArr)   
    for i = 1:1:size(points,3)
        if(i==1)
           [xunit,yunit] = circlexy(stemcenter(1,1),stemcenter(1,2),RArr(1));
         %  points(xunit,yunit,i) = 1;
        end
        if(i==size(points,3))
           [xunit,yunit] = circlexy(stemcenter(end,1),stemcenter(end,2),RArr(end));
%           points(xunit,yunit,i) = 1;
        end
    end

    
%     for i = 1:1:size(points,3)
%         valll = points(:,:,i);
%         if(or( and(sum(valll(:))>1,prev==0) , and(sum(valll(:))==0,prev>0) ))
%             [roww,coll] = ind2sub([size(points,1), size(points,2)],find(points(:,:,i)>-1)); % meshgrid([1:1:75],[1:1:75]);
% 
%             if(cnt == 0)
%                 cnt = cnt +1;
%                 [xunit1,yunit1] = ind2sub([size(points,1), size(points,2)],find(points(:,:,i)>0));
%             end
% 
%             [in,~] = inpolygon(roww(:),coll(:),xunit1, yunit1);
%             [rd,cd] = ind2sub([size(points,1), size(points,2)],find(in>0));
%             if(and(sum(valll(:))>1,prev==0) )
%                 points(5:20,5:20,i) = 1;
%             end
%         end
%         prev = sum(valll(:));
%     end
%     for i = 760:1:900
%         arr = [1,2; 1,3; 1,4; 1,5; 3,1; 3,5; 3,8; 8,2; 8,3; 8,5; 10,1; 10,5; 10,8];
%         points(arr(:,1),arr(:,2),i) = 1;
%     end
% points(1:10,1:10,900) = 1;
end

function [xunit,yunit] = circlexy(x,y,r)
    th = 0:pi/10:2*pi;
    xunit = round(r * cos(th) + x,0);
    yunit = round(r * sin(th) + y,0);
end