function [lines,lineImg] = FindParallelLines(threshPoleImg, signCoords)

GxR = [-1 0 1];
GxL = [1 0 -1];

lineRImg = conv2(GxR, single(threshPoleImg));
lineLImg = conv2(GxL, single(threshPoleImg));

lineImg = lineRImg | lineLImg;

count = 0;
startY = min(signCoords(:,1));
startX = [min(signCoords(:,2)), max(signCoords(:,2))];

lLine = lineImg(startY:end, startX(1));
rLine = lineImg(startY:end, startX(2));

lCount = sum(lLine);
rCount = sum(rLine);

window = 30;
lCol = startX(1);
rCol = startX(2);

for i =-window:window
    
    if(i < 0)
        tmpLCount = sum(lineImg(startY:end, startX(1) + i));
        
        if(tmpLCount > lCount)
            lCol = startX(1) + i;
        end
    elseif(i == 0)
        tmpLCount = sum(lineImg(startY:end, startX(1) + i));
        tmpRCount = sum(lineImg(startY:end, startX(2) + i));
        
        if(tmpRCount > rCount)
            rCol = startX(2) + i;
        end
        if(tmpLCount > lCount)
            lCol = startX(1) + i;
        end
    else
            tmpRCount = sum(lineImg(startY:end, startX(2) + i));
            
            if(tmpRCount > rCount)
                rCol = startX(2) + i;
            end
    end
end

endLY = startY + find(lineImg(startY:end, lCol), 1, 'last' ) - 1;
endRY = startY + find(lineImg(startY:end, rCol), 1, 'last' ) - 1;

lLine = [lCol, lCol; startY, endLY ];
rLine = [rCol, rCol; startY, endRY ];

num = length(find(lineImg));

[I, J] = ind2sub(size(lineImg), find(lineImg));

iter = 4;
data = [J'; I'];

%threshDist = floor((startX(2) - startX(1))/2) - 1;
threshDist = 2;

% bestInNum = 0;
% 
% finalStartPoint = [];
% finalEndPoint   = [];

startPoints = [lLine(:, 1), rLine(:, 1)];
endPoints   = [lLine(:, 2), rLine(:, 2)];

poleLines = [];

for k = 1:length(startPoints)
    
    startPoint = startPoints(:, k);
    endPoint   = endPoints(:, k);
    
    for i = 1:iter
        
        kLine = endPoint - startPoint;
        kLineNorm = kLine/norm(kLine);
        
        normVector = [-kLineNorm(2), kLineNorm(1)];
        
        distance = normVector * (data - repmat(startPoint, 1, num));
        
        inlierIdx = find(abs(distance) <= threshDist);
        inlierNum = length(inlierIdx);
        
        f_x = data(1, inlierIdx);
        f_y = data(2, inlierIdx);
        
        B = [ones(inlierNum, 1) f_y'] \ f_x';
        
        startPoint = [B(1); startPoint(2)];
        endPoint   = [B(1); endPoint(2)];
    end
    
    poleLines = [poleLines, B];
end

% y = 1:size(lineImg, 1);
% xL = poleLines(2,1)*y + poleLines(1,1);
% xR = poleLines(2,2)*y + poleLines(1,2);

lines = poleLines;

% figure
% imshow(lineImg, [])
% hold on
% plot(xL, y, 'LineWidth', 2)
% plot(xR, y, 'LineWidth', 2)
% hold off

end