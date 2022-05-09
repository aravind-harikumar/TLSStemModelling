function [c,r] = plotCircles(InputImage, ploton)
    % read image
    %InputImage= imread('Circle.jpg');
    bwInputImage = im2bw(InputImage,graythresh(InputImage));
    %figure;
    %subplot(1,2,1);
    %imshow(bwInputImage);
    % clear small aconnected components in the
    %bwInputImage = medfilt2(bwInputImage);
    bw = bwareaopen(bwInputImage,9);
    [~, threshold] = edge(bw, 'sobel');
    fudgeFactor = .5;
    im = edge(bw,'sobel', threshold * fudgeFactor);
    %subplot(1,2,2);
    if(ploton)
        hold on;
        %imshow(im);
        hold on;
    end
    % finds circles in radius range
    [c,r] = imfindcircles(im,[4,15]);
    if(ploton)
        if(length(c)>0)
            gg = c;
            gg(1) = gg(1)+1;
            gg(2) = gg(2)-0
            viscircles(gg,r);
        end
    end
end