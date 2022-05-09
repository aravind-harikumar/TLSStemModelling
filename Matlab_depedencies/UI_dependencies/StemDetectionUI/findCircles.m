function [c,r] = findCircles(InputImage,visCir)
%     close all;
    %binarization using Otsu's threshold
    bwInputImage = imbinarize(InputImage,graythresh(InputImage));
%     bwInputImage = histeq(InputImage);
    %figure;
%     subplot(1,2,1);
    
    %clear small aconnected components in the
    bw = bwareaopen(bwInputImage,9);
    [~, threshold] = edge(bw, 'sobel');
    fudgeFactor = .5;
    bw = edge(bw,'sobel', threshold * fudgeFactor);
%     subplot(1,2,2);
%     close all;
%     imshow(im);hold on;
    %finds circles in radius range
%     imshow(bw);
    [c,r] = imfindcircles(bwInputImage,[6 40]);
    if(visCir)
        hold on;
        viscircles(c,r);
        hold off;
    end
end