function [minCol, maxCol, signCoords, hVal, sVal, iVal] = checkBelowArea(jointCandidate, hsvImg)
hVal=[];
sVal=[];
iVal=[];
minCol=[];
maxCol=[];
signCoords=[];
if ~isempty(jointCandidate)
    for i = 1:size(jointCandidate,2)
        if isempty(jointCandidate{i})
            continue;
        end
        rIdx = jointCandidate{i};
        [Y,X] = ind2sub([size(hsvImg, 1) size(hsvImg, 2)], rIdx);

        winSize = 40;
        minCol = min(X) - winSize;
        maxCol = max(X) + winSize;

        if(minCol < 1)
            minCol = 1;
        end

        if(maxCol > size(hsvImg, 2))
            maxCol = size(hsvImg, 2);
        end
        
        
        wBelow = [max(Y)+5 , max(Y) + 30, min(X), max(X)];       
        
        hVals = hsvImg(wBelow(1):wBelow(2), wBelow(3):wBelow(4), 1);
        hVal = mean(mean(hVals));
        
        sVals = hsvImg(wBelow(1):wBelow(2), wBelow(3):wBelow(4), 2);
        sVal = mean(mean(sVals));
       
        iVals = hsvImg(wBelow(1):wBelow(2), wBelow(3):wBelow(4), 3);
        iVal = mean(mean(iVals));
        
        if min(Y) > 30
            wAbove = [min(Y)-30 , min(Y) - 5, min(X), max(X)];
            
            hValsA = hsvImg(wAbove(1):wAbove(2), wAbove(3):wAbove(4), 1);
            hValA = mean(mean(hValsA));
            
            sValsA = hsvImg(wAbove(1):wAbove(2), wAbove(3):wAbove(4), 2);
            sValA = mean(mean(sValsA));
            
            iValsA = hsvImg(wAbove(1):wAbove(2), wAbove(3):wAbove(4), 3);
            iValA = mean(mean(iValsA));
            
%             [hValA sValA iValA; hVal sVal iVal]
%             figure(20), imshow(hsv2rgb(hsvImg))
%             Aplot = [wAbove(3) wAbove(1) wAbove(4)-wAbove(3) wAbove(2)-wAbove(1)];
%             Bplot = [wBelow(3) wBelow(1) wBelow(4)-wBelow(3) wBelow(2)-wBelow(1)];
%                     imrect(gca,Aplot);
%                     imrect(gca,Bplot);
            if (hVal > 0.4 && sVal < 0.40 && iVal < 0.40) && (abs(hVal-hValA) > 0.1 || (abs(sVal-sValA) > 0.1 && abs(iVal-iValA) > 0.1)) % Ant iVal < 0.30
%             if (sVal < 0.40 && iVal < 0.40) && (abs(hVal-hValA) > 0.1 || (abs(sVal-sValA) > 0.1 || abs(iVal-iValA) > 0.1)) % Ant iVal < 0.30
                signCoords = [Y, X];
                break;
            end
        else 
%             figure(20), imshow(hsv2rgb(hsvImg))
%             Bplot = [wBelow(3) wBelow(1) wBelow(4)-wBelow(3) wBelow(2)-wBelow(1)];
%             imrect(gca,Bplot);
%             [hVal sVal iVal]
            if (hVal > 0.4 && sVal < 0.40 && iVal < 0.40) % Ant iVal < 0.30
%             if (sVal < 0.40 && iVal < 0.40)
                signCoords = [Y, X];
                    break;
            end
        end 
    end  
end