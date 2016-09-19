function [output , match] = checkSplitTrunk(CC2, thr, mid)

%treeList = List of connected components identified as tree cluster
%thr: Threshold number to select the amount of connected components neighbours (+30 pixels) wanted to check.
%mid: horitzontal middle of the picture.
%Sergi Salgueiro
treeList = CC2.PixelIdxList;
possComp=[];
match=[];
    for i=1:length(treeList)
        [Y,X]=ind2sub(CC2.ImageSize,treeList{i});

        if(size(treeList{i}, 1) > 30) && mean(Y) <= mid
            cornersPossComp=ApproxRectangle([X,Y]);
            centersPossComp=[(cornersPossComp(1,1)+cornersPossComp(3,1))/2 (cornersPossComp(1,2)+cornersPossComp(2,2))/2];
            possComp=[possComp,i];

            if length(possComp)-thr > 1
                k=length(possComp)-thr;
            else
                k=1;
            end

            if i > 1
            for j=k:length(possComp)-1
                if ~isempty(treeList{possComp(j)})
                    [Y2,X2]=ind2sub(CC2.ImageSize,treeList{possComp(j)});
                    cornersEv=ApproxRectangle([X2,Y2]);
                    centersEv=[(cornersEv(1,1)+cornersEv(3,1))/2 (cornersEv(1,2)+cornersEv(2,2))/2];            
                    % Comprova que els troncs estan a la mateixa vertical
                    possCompInEv = cornersEv(1,1) <= centersPossComp(1) && centersPossComp(1) <= cornersEv(3,1);
                    EvInPossComp = cornersPossComp(1,1) <= centersEv(1) && centersEv(1) <= cornersPossComp(3,1);
                    % Comprova que els elements estan un sobre l'altre
                    yAvg=mean(Y2)- mean(Y); % positiu=> Ev sota possComp. Negatiu=>possComp sota Ev
                    if yAvg >= 0
                        vertCon = max(Y) < min(Y2);
                    else
                        vertCon = max(Y2) < min(Y);
                    end

                    if (possCompInEv || EvInPossComp) && vertCon
                        match = [match;possComp(j) i];
                        treeList{i}=sort([treeList{possComp(j)};treeList{i}]);
                        treeList{possComp(j)}=[];
                    end
                end
            end
            end
        end

    end

output = treeList;
end