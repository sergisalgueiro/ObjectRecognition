function output = FindLightPost(farDl)
    output = [];
    farDl = bwareaopen(farDl,30);
    CC = bwconncomp(farDl,4);
    if CC.NumObjects > 1
        [LPList , ~] = checkSplitTrunk(CC, 1, ceil(size(farDl,1)*2/3));
    else 
        LPList = CC.PixelIdxList;
    end
    
    for i = 1:length(LPList)
        if isempty(LPList{i})
            continue;
        end
        [Y,X]=ind2sub(CC.ImageSize,LPList{i});
        yEnds=[min(Y) max(Y)];
        xEnds=[min(X) max(X)];
        yEndsDiff=yEnds(2)-yEnds(1);
        xEndsDiff=xEnds(2)-xEnds(1);
        if xEndsDiff < 15 && xEndsDiff > 2 && (yEndsDiff/xEndsDiff) > 9 && yEndsDiff > size(farDl,1)/2
            coordsLightPost = [xEnds(1) , yEnds(1), xEndsDiff, yEndsDiff];
            output = coordsLightPost;
            break;
        end
    end
    
end
