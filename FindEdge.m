function edgeRoute = FindEdge(roadImg)

IND = find(roadImg);

[r, c] = size(roadImg);
%Extreu X,Y dels pixels que formen carretera.
[I,J] = ind2sub(size(roadImg),IND);

markedPix = 0;
startPix = [I(1) J(1)];
%defineix quins son els pixels de sota i de sobre
if(startPix(1) + 1 > r)
    downPix = 0;
    upPix = roadImg(startPix(1) - 1, startPix(2));
elseif(startPix(1) - 1 < 1)
    upPix = 0;
    downPix = roadImg(startPix(1) + 1, startPix(2));
else
    downPix = roadImg(startPix(1) + 1, startPix(2));
    upPix = roadImg(startPix(1) - 1, startPix(2));
end
%defineix esquerra i dreta
if(startPix(2) + 1 > c)
    rightPix = 0;
    leftPix = roadImg(startPix(1), startPix(2) - 1);
elseif(startPix(2) - 1 < 1)
    leftPix = 0;
    rightPix = roadImg(startPix(1), startPix(2) + 1);
else
    leftPix = roadImg(startPix(1), startPix(2) - 1);
    rightPix = roadImg(startPix(1), startPix(2) + 1);
end

endPix = [0 0];
curPix = startPix;
count = 0;
lastPos = 3;
first = 0;
route = [];

while(markedPix ~= 1)
    count = count + 1;
    
    if(curPix(1) + 1 > r)
        down = curPix;
        downPix = 0;
        
        up = [curPix(1) - 1, curPix(2)];
        upPix = roadImg(up(1), up(2));
    elseif(curPix(1) - 1 < 1)
        up = curPix;
        upPix = 0;
        
        down = [curPix(1) + 1, curPix(2)];
        downPix = roadImg(down(1), down(2));
    else
        down = [curPix(1) + 1, curPix(2)];
        downPix = roadImg(down(1), down(2));
        
        up = [curPix(1) - 1, curPix(2)];
        upPix = roadImg(up(1), up(2));
    end
    
    if(curPix(2) + 1 > c)
        right = curPix;
        rightPix = 0;
        
        left = [curPix(1), curPix(2) - 1];
        leftPix = roadImg(left(1), left(2));
    elseif(curPix(2) - 1 < 1)
        left = curPix;
        leftPix = 0;
        
        right = [curPix(1), curPix(2) + 1];
        rightPix = roadImg(right(1), right(2));
    else
        left = [curPix(1), curPix(2) - 1];
        leftPix = roadImg(left(1), left(2));
        
        right = [curPix(1), curPix(2) + 1];
        rightPix = roadImg(right(1), right(2));
    end
    
    if(all(curPix == endPix) && first == 1)
        break;
    end
    
    if(first == 0)
        if(leftPix == 1)
            curPix = left;
            lastPos = 2;
        else
            if(upPix == 1)
                endPix = curPix;
                curPix = up;
                lastPos = 3;
                
                route = [route; endPix, lastPos];
            elseif(rightPix == 1)
                endPix = curPix;
                curPix = right;
                lastPos = 4;
                
                route = [route; endPix, lastPos];
            elseif(downPix == 1)
                endPix = curPix;
                curPix = down;
                lastPos = 1;
                
                route = [route; endPix, lastPos];
            end
            
            route = [route; curPix, lastPos];
            
            first = 1;
        end
    else
        if(lastPos == 1)
            if(rightPix == 1)
                curPix = right;
                lastPos = 4;
            elseif(downPix == 1)
                curPix = down;
                lastPos = 1;
            elseif(leftPix == 1)
                curPix = left;
                lastPos = 2;
            elseif(upPix == 1)
                curPix = up;
                lastPos = 3;
            end
        elseif(lastPos == 2)
            if(downPix == 1)
                curPix = down;
                lastPos = 1;
            elseif(leftPix == 1)
                curPix = left;
                lastPos = 2;
            elseif(upPix == 1)
                curPix = up;
                lastPos = 3;
            elseif(rightPix == 1)
                curPix = right;
                lastPos = 4;
            end
        elseif(lastPos == 3)
            if(leftPix == 1)
                curPix = left;
                lastPos = 2;
            elseif(upPix == 1)
                curPix = up;
                lastPos = 3;
            elseif(rightPix == 1)
                curPix = right;
                lastPos = 4;
            elseif(downPix == 1)
                curPix = down;
                lastPos = 1;
            end
        elseif(lastPos == 4)
            if(upPix == 1)
                curPix = up;
                lastPos = 3;
            elseif(rightPix == 1)
                curPix = right;
                lastPos = 4;
            elseif(downPix == 1)
                curPix = down;
                lastPos = 1;
            elseif(leftPix == 1)
                curPix = left;
                lastPos = 2;
            end
        end
        
        route = [route; curPix, lastPos];
    end
end

edgeRoute = route;

end