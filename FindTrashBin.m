function [tbImg,box] = FindTrashBin(tbImg)
    % OPEN CLOSE ultim moment
%     se = strel('rectangle', [3 3]);
%     tbEr = imerode(tbImg, se);
%     tbImg = imdilate(tbEr, se);
    %
    tbImg = bwareaopen(tbImg,40);
    tbImg = imfill(tbImg,'holes');
    stats = regionprops(tbImg,'BoundingBox','Area');
    box=[];
    for i = 1 : length(stats)
        dY=stats(i).BoundingBox(4);
        dX=stats(i).BoundingBox(3);
        areaRectangle=dX*dY;
        areaRate=stats(i).Area/areaRectangle;
        yxRate=dY/dX;
        if areaRate >= 0.60  && (yxRate >= 1.5 && yxRate <= 2.3)
           box=stats(i).BoundingBox;
           break;
        end
    end
end