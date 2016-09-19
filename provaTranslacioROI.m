piCamRange = 53.5;
if  strcmp('red',odoCase) || strcmp('red2',odoCase)
ROI=[tbBox(1)-tbBox(3)*2 tbBox(2)-(tbBox(4)/2) round(tbBox(3)*5) round(tbBox(4)*2)];
else
 ROI=[tbBox(1)-tbBox(3) tbBox(2)-(tbBox(4)/2) round(tbBox(3)*3) round(tbBox(4)*2)];
end

object=dImg(tbBox(2):tbBox(2)+tbBox(4),tbBox(1):tbBox(1)+tbBox(3),3);
% Provant el match frame to frame amb SIFT
% grImg = rgb2gray(sImg);
% grObject=grImg(tbBox(2)-10:tbBox(2)+tbBox(4)+10,tbBox(1)-10:tbBox(1)+tbBox(3)+10);
%
figure
imshow(dImg)
imrect(gca,tbBox);
%%
S2 = num2str(nIt);
nImgName=strcat(S1,S2,S3);
tImg = imread(nImgName);
Img = tImg(:, 1:end-4, :);
scale = 1/5;
sImg = imresize(Img, scale);
hImg = imresize(Img, 1/2);
if  strcmp('red',odoCase) || strcmp('red2',odoCase)
    dImg = hImg;
else
    dImg = sImg;
end
%%
nX = (size(dImg,2)/2) - (floor(rad2deg(newCamAngle)*(size(dImg,2)/2)/(piCamRange/2)+(ROI(3)/2)));
if nX < 1
    nX=1;
end
% figure;imshow(dImg)
% imrect(gca,[nX ROI(2:4)]);
% % imrect(gca,[nX tbBox(2:4)]);

%Correction for ROI that go beyond the image borders
if nX+ROI(3) > size(dImg,2)
    x_ax = size(dImg,2);
else
    x_ax = nX+ROI(3);
end
if ROI(2) < 1
    ROI(4) = ROI(4) + ROI(2);
    ROI(2) = 1;
end

vImg=dImg(ROI(2):ROI(2)+ROI(4),nX:x_ax,3);
if size(vImg,1) > size(object,1) && size(vImg,2) > size(object,2)
    c = normxcorr2(object,vImg);
    figure, surf(c), shading flat
    xlabel('ROI X dimension (px)')
    ylabel('ROI Y dimension (px)')
    zlabel('NCC result')

    % Find peak in cross-correlation.
    [ypeak, xpeak] = find(c==max(c(:)));
   
    % Account for the padding that normxcorr2 adds.
    yoffSet = ypeak-size(object,1);
    xoffSet = xpeak-size(object,2);

    % Display matched area.
%     hFig = figure;
%     hAx  = axes;
%     imshow(vImg,'Parent', hAx);
%     imrect(hAx, [xoffSet, yoffSet, size(object,2), size(object,1)]);
    boxNCC = [xoffSet+nX, yoffSet+ROI(2), size(object,2), size(object,1)];
%     imrect(hAx, ROI);
end