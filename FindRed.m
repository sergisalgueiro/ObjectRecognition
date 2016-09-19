function [redimg, sign, coords, rsBox, minCol, poleLines, lineImg] = FindRed(sImg,hsvImg,mode)
sign=[]; coords=[]; rsBox=[]; minCol=[]; poleLines=[]; lineImg=[];
highThresh = 0.05;
lowThresh = 0.9;    
threshImg = (hsvImg(:,:,1) > lowThresh) | (hsvImg(:,:,1) < highThresh);

rImg = sImg(:,:,1);
gImg = sImg(:,:,2);
bImg = sImg(:,:,3);

se = strel('square', 2);
threshImg2 = threshImg;
threshImg2 = imopen(threshImg2, se);
threshImg2 = imclose(threshImg2, se);

redimg=threshImg2;

if(isempty(find(threshImg2, 1)))
    sign = [];
    coords = [];
    rsBox=[];
    poleLines = [];
    lineImg = [];    
    return;
end
%%
switch (mode)
    case 'sign'
         % METODE RGB
        thrRG = 15; thrRB = 15; maxDif = [1 1]; meanDif = [0 0];
        while meanDif(1) / maxDif(1) < 0.5 && meanDif(2) / maxDif(2) < 0.5
            Cond = rImg - gImg > thrRG & rImg-bImg > thrRB;
            if nnz(Cond) == 0;
                sign = [];
                coords = [];
                rsBox=[];
                poleLines = [];
                lineImg = [];
                minCol = [];
                return;
            end
            threshRGBImg = Cond;
            %Si un pixel vermell te 0 de valor, despres th#V no tindran la mateixa mida
            rImg(rImg == 0) = 1;
            gImg(gImg == 0) = 1;
            bImg(bImg == 0) = 1;
            %
            thr=uint8(threshRGBImg) .* rImg;
            thg=uint8(threshRGBImg) .* gImg;
            thb=uint8(threshRGBImg) .* bImg;
            th(:,:,1) = thr;
            th(:,:,2) = thg;
            th(:,:,3) = thb;
            thrV = thr(thr ~= 0);
            thgV = thg(thg ~= 0);
            thbV = thb(thb ~= 0);
            diffRG= thrV - thgV;
            diffRB= thrV - thbV;
            maxDif = [max(diffRG) max(diffRB)];
            meanDif = [mean(diffRG) mean(diffRB)];
            if (meanDif(1) / maxDif(1) < 0.6)
                thrRG = thrRG+1;
            end
            if (meanDif(2) / maxDif(2) < 0.6)
                thrRB = thrRB+1;
            end    
        end

        mergeImg = threshImg & threshRGBImg;
        redimg = mergeImg;
    % Filter shapes
    ccList = bwconncomp(mergeImg);
%     ccList = bwconncomp(threshImg2);
    signCandidate=[];
    largeCandidate=[];
    verticalCutOff=[31 190];
    for i = 1:ccList.NumObjects
        if length(ccList.PixelIdxList{i}) < 800 && length(ccList.PixelIdxList{i}) > 40
            [Y,X] = ind2sub(size(threshImg2), ccList.PixelIdxList{i});
            hwRatio=(max(Y)-min(Y))/(max(X)-min(X));
             if  (hwRatio <= 3/2 && hwRatio >= 0.63) && (mean(Y) < verticalCutOff(2))
                 signCandidate = [signCandidate, ccList.PixelIdxList(i)];
             elseif (hwRatio < 5/2 && hwRatio > 3/2) && (mean(Y) < verticalCutOff(2))
                 largeCandidate = [largeCandidate, ccList.PixelIdxList(i)];
             end
        end
    end

    if ~isempty(signCandidate) && length(signCandidate) > 1
        [jointCandidate, ~] = checkSplitRed(signCandidate, 1, ccList.ImageSize);
    else
        jointCandidate = signCandidate;
    end
    if isempty(jointCandidate) && isempty(largeCandidate)
        sign = [];
        coords = [];
        rsBox=[];
        poleLines = [];
        lineImg = [];
        minCol = [];
        return;
    end
    % Area below

            [minCol, maxCol, signCoords, hVal, sVal, iVal] = checkBelowArea(jointCandidate, hsvImg);
            if isempty(signCoords)
               [minCol, maxCol, signCoords, hVal, sVal, iVal] = checkBelowArea(largeCandidate, hsvImg);
            end   

            if isempty(signCoords)%%
                sign = [];
                coords = [];
                poleLines = [];
                lineImg = [];
                rsBox=[];
                return;
            end

            poleImg = hsvImg(:, minCol:maxCol, :);
    %         minCol
            hThresh = 0.2;
            sThresh = 0.15;
            iThresh = 0.1;

            threshPoleImg = ((poleImg(:,:,1) < hVal + hThresh) & (poleImg(:,:,1) >= hVal - hThresh)) & ...
                            ((poleImg(:,:,2) < sVal + sThresh) & (poleImg(:,:,2) >= sVal - sThresh)) & ...
                            ((poleImg(:,:,3) < iVal + iThresh) & (poleImg(:,:,3) >= iVal - iThresh));

            signCoords = [signCoords(:, 1), signCoords(:, 2) - minCol];
            threshPoleImg(signCoords(:,1), signCoords(:, 2)) = 1;

            sign = threshPoleImg;
            coords = signCoords;

    % Find paralel lines
        if(~isempty(threshPoleImg))
            [poleLines,lineImg] = FindParallelLines(threshPoleImg, signCoords);
        end
    % Threshold and Box and compare lines and Box
    se = strel('square', 1);
    threshPoleImg2= bwareaopen(threshPoleImg,50);
    threshPoleImg2 = imopen(threshPoleImg2, se);
    threshPoleImg2 = imclose(threshPoleImg2, se);
    ccTPI = bwconncomp(threshPoleImg2,8);
    [~, match] = checkSplitRed(ccTPI.PixelIdxList, 1, ccTPI.ImageSize);
    msc = mean(signCoords);
    stats = regionprops(threshPoleImg2,'BoundingBox');
    if ~isempty(match)
        xma=0;
        xmi=size(threshPoleImg2,2);
        yma=0;
        ymi=size(threshPoleImg2,1);
        for i=1:length(match)
            if stats(i).BoundingBox(1) < xmi
                xmi = stats(i).BoundingBox(1);
            end
            if stats(i).BoundingBox(2) < ymi
                ymi = stats(i).BoundingBox(2);
            end
            if stats(i).BoundingBox(1) + stats(i).BoundingBox(3) > xma
                xma = stats(i).BoundingBox(1) + stats(i).BoundingBox(3);
            end
            if stats(i).BoundingBox(2) + stats(i).BoundingBox(4) > yma
                yma = stats(i).BoundingBox(2) + stats(i).BoundingBox(4);
            end
            if i~= length(match)
                stats(i).BoundingBox = [];
            else
                stats(i).BoundingBox = [xmi ymi xma-xmi yma-ymi];
            end
        end
    end
    
    for i = 1:length(stats)
        if ~isempty(stats(i).BoundingBox)
        xCondition = stats(i).BoundingBox(1) < msc(2) &&  stats(i).BoundingBox(1)+stats(i).BoundingBox(3) > msc(2);
        yCondition = stats(i).BoundingBox(2) < msc(1) &&  stats(i).BoundingBox(2)+stats(i).BoundingBox(4) > msc(1);
        if xCondition && yCondition
            if (stats(i).BoundingBox(4)/stats(i).BoundingBox(3)) > 6
                rsBox = stats(i).BoundingBox;
                rsBox(1) = rsBox(1) + minCol-1;
            elseif (stats(i).BoundingBox(4)/stats(i).BoundingBox(3)) > 3
                rsBox = stats(i).BoundingBox;
                rsBox(1) = rsBox(1) + minCol-1;
                linear = 0;
                switch linear
                    case 1
                    interval=rsBox(2):rsBox(2)+rsBox(4);
                    xLeft = poleLines(2,1)*interval + poleLines(1,1);
                    xRight = poleLines(2,2)*interval + poleLines(1,2);

                    leftEdge = min(xLeft) + minCol -1;
                    rightEdge = max(xRight) + minCol -1;
                    if abs(leftEdge - rsBox(1)) > 5 && abs(rightEdge - (rsBox(1)+rsBox(3))) < 5
                        rsBox(3) = ceil(rsBox(1)+rsBox(3) - leftEdge);
                        rsBox(1) = ceil(leftEdge); 
                    elseif abs(leftEdge - rsBox(1)) < 5 && abs(rightEdge - (rsBox(1)+rsBox(3))) > 5
                        rsBox(3) = ceil(rightEdge - rsBox(1));
                        rsBox(1) = round(rsBox(1)); 
                    else
                        rsBox=[];
                        sign = [];
                        coords = [];
                        poleLines = [];
                        lineImg = [];
                    end
                    case 0
                    %APPLYING DISTANCE THRESHOLD
                    leftEdge = min(signCoords(:,2))+ minCol -1;
                    rightEdge = max(signCoords(:,2))+ minCol -1;
                    upperEdge = min(signCoords(:,1));
                    if   abs(leftEdge - rsBox(1)) > 5
                        rsBox(3) = ceil(rsBox(1)+rsBox(3) - leftEdge);
                        rsBox(1)=leftEdge;
                    end
                    if abs(rightEdge - (rsBox(1)+rsBox(3))) > 5
                        rsBox(3) = rightEdge - rsBox(1)+3;
                    end
                    if abs(upperEdge - rsBox(2)) > 5
                        rsBox(3) = upperEdge - 1;
                    end
                 end
            else
                rsBox=[];
                sign = [];
                coords = [];
                poleLines = [];
                lineImg = [];
            end
        end
        end
    end
%%
    case 'mid'
        sign = [];
        coords = [];
        rsBox=[];
        minCol=[];
 
        % METODE RGB
        thrRG = 15; thrRB = 15; maxDif = [1 1]; meanDif = [0 0];
        while meanDif(1) / maxDif(1) < 0.5 && meanDif(2) / maxDif(2) < 0.5
            Cond = rImg - gImg > thrRG & rImg-bImg > thrRB;
            if nnz(Cond) == 0;
                sign = [];
                coords = [];
                rsBox=[];
                poleLines = [];
                lineImg = [];
                minCol = [];
                return;
            end
            threshRGBImg = Cond;
            %Si un pixel vermell te 0 de valor, despres th#V no tindran la mateixa mida
            rImg(rImg == 0) = 1;
            gImg(gImg == 0) = 1;
            bImg(bImg == 0) = 1;
            %
            thr=uint8(threshRGBImg) .* rImg;
            thg=uint8(threshRGBImg) .* gImg;
            thb=uint8(threshRGBImg) .* bImg;
            th(:,:,1) = thr;
            th(:,:,2) = thg;
            th(:,:,3) = thb;
            thrV = thr(thr ~= 0);
            thgV = thg(thg ~= 0);
            thbV = thb(thb ~= 0);
            diffRG= thrV - thgV;
            diffRB= thrV - thbV;
            maxDif = [max(diffRG) max(diffRB)];
            meanDif = [mean(diffRG) mean(diffRB)];
            if (meanDif(1) / maxDif(1) < 0.6)
                thrRG = thrRG+1;
            end
            if (meanDif(2) / maxDif(2) < 0.6)
                thrRB = thrRB+1;
            end    
        end

        mergeImg = threshImg & threshRGBImg;
        redimg = mergeImg;
        % Identification/Segmentation
        mergeImg = bwareaopen(mergeImg,50);
        ccList = bwconncomp(mergeImg);
        midSignNum = [];
        for i=1:ccList.NumObjects
            [Y,X] = ind2sub(size(threshImg2), ccList.PixelIdxList{i});
            yxRatio = (max(Y)-min(Y))/(max(X)-min(X));
            areaRatio = length(ccList.PixelIdxList{i})/((max(Y)-min(Y))*(max(X)-min(X)));
            if yxRatio > 0.5 && yxRatio < 0.75 && size(ccList.PixelIdxList{i},1) > 100 && areaRatio > 0.3
                tempBox = [min(X) min(Y) max(X)-min(X) max(Y)-min(Y)];
                midSignNum = i;
                break;
            end 
        end
        if ~isempty(midSignNum)
            h_Img = hsvImg(:,:,3);
            yStart = min(Y)-10;
            if yStart <1
                yStart = 1;
            end
            xStart = min(X)-10;
            if xStart <1
                xStart = 1;
            end 
            xEnd = max(X)+10;
            if xEnd > size(h_Img,2)
                xEnd = size(h_Img,2);
            end            
            midSign = h_Img(yStart:end,xStart:xEnd,:);
            center = [tempBox(1)+tempBox(3)/2 tempBox(2)+tempBox(4)/2];
            postBelow = edge(midSign(max(Y)-min(Y)+10:floor(max(Y)-min(Y)+10+(max(Y)-min(Y))*1.5),:),'Canny');
            lineImg = postBelow;
            % try hough
            [H,T,R] = hough(postBelow);
            P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
            % adapting vertical lines
            P1=P(:,2);
            P1=P1(P1>80 & P1<100);
            P=P(P(:,2) > 80 & P(:,2) < 100);
            P=[P P1];
            lines = houghlines(postBelow,T,R,P,'FillGap',30,'MinLength',2);
            poleLines = lines;
            
            if ~isempty(P)
                xcenter = round(center(1)-xStart);
                linesL = [];
                linesR = [];

                for i=1:length(lines)
                    locCond = [lines(i).point1(1) < xcenter && lines(i).point2(1) < xcenter lines(i).point1(1) > xcenter && lines(i).point2(1) > xcenter]; 
                    thCond = lines(i).theta < 10 && lines(i).theta > -10;
                    longCond = lines(i).point2(2) - lines(i).point1(2) >= 0.65*size(postBelow,1);
                    if locCond(1) && thCond && longCond
                        linesL = [linesL; lines(i).point1 lines(i).point2 lines(i).theta];
                    elseif locCond(2) && thCond && longCond
                        linesR = [linesR; lines(i).point1 lines(i).point2 lines(i).theta];
                    end
                end

                if ~isempty(linesR) && ~isempty(linesL) && (size(linesL,1)+size(linesR,1)) >= 3
                    %Most probable detection
                    rsBox = tempBox;
                else
                    rsBox = [];
                    lineImg = [];
                end        
            end
        else 
            poleLines = [];
            rsBox = [];
            lineImg = [];
        end
        
   %%
    case 'big'
        sign = [];
        coords = [];
        rsBox=[];
        minCol=[];
        thrRG = 15; thrRB = 15; maxDif = [1 1]; meanDif = [0 0];
        while meanDif(1) / maxDif(1) < 0.5 && meanDif(2) / maxDif(2) < 0.5
            Cond = rImg - gImg > thrRG & rImg-bImg > thrRB;
            if nnz(Cond) == 0;
                sign = [];
                coords = [];
                rsBox=[];
                poleLines = [];
                lineImg = [];
                minCol = [];
                return;
            end
            threshRGBImg = Cond;
            %Si un pixel vermell te 0 de valor, despres th#V no tindran la mateixa mida
            rImg(rImg == 0) = 1;
            gImg(gImg == 0) = 1;
            bImg(bImg == 0) = 1;
            %
            thr=uint8(threshRGBImg) .* rImg;
            thg=uint8(threshRGBImg) .* gImg;
            thb=uint8(threshRGBImg) .* bImg;
            th(:,:,1) = thr;
            th(:,:,2) = thg;
            th(:,:,3) = thb;
            thrV = thr(thr ~= 0);
            thgV = thg(thg ~= 0);
            thbV = thb(thb ~= 0);
            diffRG= thrV - thgV;
            diffRB= thrV - thbV;
            maxDif = [max(diffRG) max(diffRB)];
            meanDif = [mean(diffRG) mean(diffRB)];
            if (meanDif(1) / maxDif(1) < 0.6)
                thrRG = thrRG+1;
            end
            if (meanDif(2) / maxDif(2) < 0.6)
                thrRB = thrRB+1;
            end    
        end
        mergeImg = threshImg & threshRGBImg;
        redimg = mergeImg;
        % Identification/Segmentation
        mergeImg = bwareaopen(mergeImg,50);
        
        se = strel('square', 5);
        promergeImg = imopen(mergeImg, se);
        promergeImg = imclose(promergeImg, se);

        ccList = bwconncomp(promergeImg);
        for i=1:ccList.NumObjects
            [Y,X] = ind2sub(size(promergeImg), ccList.PixelIdxList{i});
            yxRatio = (max(Y)-min(Y))/(max(X)-min(X));
            if yxRatio > 1.1 && yxRatio < 1.3 && size(ccList.PixelIdxList{i},1) > 750
                img = zeros(size(threshImg));
                img(ccList.PixelIdxList{i}) = 1;
                [~,detect]=FindRectangle(img);
                if detect
                    rsBox = [min(X) min(Y) max(X)-min(X) max(Y)-min(Y)];
                    break;
                end
            end 
        end
end
end