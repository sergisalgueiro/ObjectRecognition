function [output , match] = checkSplitRed(CC2, thr, imsize)

%treeList = List of connected components identified as tree cluster
%thr: Threshold number to select the amount of connected components neighbours (+30 pixels) wanted to check.
%mid: horitzontal middle of the picture.
%Sergi Salgueiro
redList = CC2;
possComp=[];
match=[];
    for i=1:length(redList)
        [Y,X]=ind2sub(imsize,redList{i});

%             cornersPossComp=ApproxRectangle([X,Y]);
%             centersPossComp=[(cornersPossComp(1,1)+cornersPossComp(3,1))/2 (cornersPossComp(1,2)+cornersPossComp(2,2))/2];
            possComp=[possComp,i];

            if length(possComp)-thr > 1
                k=length(possComp)-thr;
            else
                k=1;
            end

            if i > 1
            for j=k:length(possComp)-1
                if ~isempty(redList{possComp(j)})
                    [Y2,X2]=ind2sub(imsize,redList{possComp(j)});
%                     cornersEv=ApproxRectangle([X2,Y2]);
%                     centersEv=[(cornersEv(1,1)+cornersEv(3,1))/2 (cornersEv(1,2)+cornersEv(2,2))/2];            
                    % Comprova que els troncs estan a la mateixa vertical
                    possCompInEv =  min(X) <= min(X2) && max(X) >= max(X2);%cornersEv(1,1) <= centersPossComp(1) && centersPossComp(1) <= cornersEv(3,1);
                    EvInPossComp = min(X2) <= min(X) && max(X2) >= max(X);%cornersPossComp(1,1) <= centersEv(1) && centersEv(1) <= cornersPossComp(3,1);
                    % Comprova que els elements estan un sobre l'altre
                    yAvg=mean(Y2)- mean(Y); % positiu=> Ev sota possComp. Negatiu=>possComp sota Ev
                    if yAvg >= 0
                        vertCon = 10 > min(Y2) - max(Y);
                    else
                        vertCon = 10 > min(Y) - max(Y2);
                    end

                    if (possCompInEv || EvInPossComp) && vertCon
                        match = [match;possComp(j) i];
                        redList{i}=sort([redList{possComp(j)};redList{i}]);
                        redList{possComp(j)}=[];
                    end
                end
            end
            end
    end

output = redList;
end