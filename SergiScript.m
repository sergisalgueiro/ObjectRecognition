% clc
close all
clear all
% first = true;
plots=true;
road=1;
tree=1;
redSign=0;
midSign = 0;
bigSign = 0;
trashbin=0;
lightpost=0;
%% Import images
tic
tImg = imread('l10.bmp');
Img = tImg(:, 1:end-4, :);
scale = 1/5;
sImg = imresize(Img, scale);
hImg = imresize(Img, 1/2);

if trashbin
    rImg = sImg(:,:,1);
    gImg = sImg(:,:,2);
    bImg = sImg(:,:,3);
end
%% HSV Tranformation
hsvImg = rgb2hsv(sImg);
hsvImg2 = rgb2hsv(hImg);
mid = size(hsvImg, 1)/2;
lenLayer = size(hsvImg,1)*size(hsvImg,2);
%cada H S V els ordena per columnes
H = hsvImg(:,:,1);
S = hsvImg(:,:,2);
V = hsvImg(:,:,3);
rawInput=[H(:) S(:) V(:)];
%% 4D Histogram
%Input �s la Imatge i l'interval
[hist3D, space, pixCor] = Hist3D(hsvImg, 25);

%% Adaptive Threshold
mH = mean(H(:));
mS = mean(S(:));
mV = mean(V(:));

w1=0.6; w2=0.4;
Hma = 0.3; Hmi = 0.16;
Sma = 0.28; Smi = 0.12;
if mH < Hmi && mS > Sma
    clustThresh = 10;
elseif mH > Hma && mS < Smi
    clustThresh = 50;
elseif mS < Smi
    clustThresh = 50;
elseif mS > Sma
    clustThresh = 20;
else
    eqH = w1 * (1-(Hma - mH)/(Hma-Hmi));
    eqS = w2 * (Sma - mS)/(Sma-Smi);
    clustThresh = 20 + round(30*(eqH+eqS));
end
%% Seeding

disp('clustThresh')
disp(clustThresh)
% clustThresh = 50;
clustCoords = [];
clustVal = clustThresh + 1;
oldClustVal = 0;
count = 1;

peaks3D = FindProminentPeak4D(hist3D, 1, 2);

while(true)
    
    if(count == 1)
        %com ones(), pero el valor �s l�gic al workspace
        peakMask = true(length(peaks3D) , 1);
        
    else           
        peakMask = (peaks3D(:, 1) > (coords(1, 1) + 1) | peaks3D(:, 1) < (coords(1, 1) - 1) )| ...
                   (peaks3D(:, 2) > (coords(1, 2) + 1) | peaks3D(:, 2) < (coords(1, 2) - 1) )| ...
                   (peaks3D(:, 3) > (coords(1, 3) + 1) | peaks3D(:, 3) < (coords(1, 3) - 1));

        peakMask = peakMask & oldMask;
    end
    
    newPeaks = peaks3D(peakMask, :);
    
    clustVal = max(newPeaks(:,4));

    if(clustVal < clustThresh)
        break;
    end
    %Enganxa les coordenades del colega amb valor m�s alt
    coords = newPeaks(newPeaks(:,4) == max(newPeaks(:,4)), 1:3);
    
    clustCoords = [clustCoords; coords];
    
    oldMask = peakMask;
    
    oldClustVal = clustVal;
    
    count = count + 1;
end
%% K-Means clustering

seed = clustCoords;
cluster = size(seed, 1);
[idx, cMat, sumd] = kmeans(rawInput, cluster, 'Start', seed./25, 'emptyaction','singleton');

%% Plotting of clusters

[hx, hy] = ind2sub([size(hsvImg, 1) size(hsvImg, 2)], 1:lenLayer);
mMat = zeros(cluster, 2);
clustImgs = cell(1, cluster);
rgbMat = zeros(cluster, 3);

for cClust = 1:cluster
    clustImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
    %Extreu X i Y de cada punt/p�xel del cluster
    cXVal = hx(idx == cClust);
    cYVal = hy(idx == cClust);
    %assigna 1 al cluster 0 a la resta, per la binaritzaci� de la gr�fica
    
    r=0; g=0; b=0;
    for i = 1:length(cXVal)
        clustImg(cXVal(i), cYVal(i)) = 1;
        if trashbin
            r=r+double(rImg(cXVal(i), cYVal(i)));
            g=g+double(gImg(cXVal(i), cYVal(i)));
            b=b+double(bImg(cXVal(i), cYVal(i)));
        end
    end
    if trashbin
        rgbMat(cClust,:) = [r/length(cXVal) g/length(cXVal) b/length(cXVal)];
    end
    %Guarda el cl�ster a la seva corresponent posici� de la cel�la    
    clustImgs{cClust} = clustImg;
    %Coordenada mitja dels pixels en la imatge
    mMat(cClust, :) = [round(sum(cXVal)/length(cXVal)), round(sum(cYVal)/length(cYVal))];
    
end

%% Classification
 
clusters = 1:cluster;

%Associa el cel al V m�s alt de tots
maxSkyClust = find(cMat(:, 3) == max(cMat(:, 3)));

skyClust = [];
%Si la saturaci� esta per sota de 0.1 i la mitjana de pixels per sobre la
%meitat de la foto, l'assigna i l'esborra.
if(cMat(maxSkyClust, 2) < 0.1 && mMat(maxSkyClust, 1) < mid)
    skyClust = maxSkyClust;
%     strcat('The sky cluster is: #', num2str(maxSkyClust))
    clusters(skyClust) = [];
end
%Calcula ratios
hvRatio = cMat(:, 1) ./ cMat(:, 3);
svRatio = cMat(:, 2) ./ cMat(:, 3);
hsRatio = cMat(:, 1) ./ cMat(:, 2);

roadClustH = [];
treeClustH = [];
trashBin=[];
lightPost=[];
maxRoadSV = zeros(1,2);
msg = '';

for i = clusters
%   usa un ratio i un m�x de value per la carretera
     if(svRatio(i) < 0.5 && cMat(i, 3) < 0.8)
            if(cMat(i, 3) > 0.40 && mMat(i, 1) >= mid)
                %aqui es queda amb el ratio mes alt dels clusters
                if(svRatio(i) > maxRoadSV(1))
                    maxRoadSV(1) = svRatio(i);
                    maxRoadSV(2) = i;
                end
                %es queda amb els que passen el tall com a molt probables
                roadClustH = [roadClustH, i];
                msg = sprintf('%s \n Cluster #%d high road prob.', msg, i);
            end
     end
    %Class paperera 
    noTree = true; 
    rCond = rgbMat(i,1) < 70;
    gCond = rgbMat(i,2) < 70;
    bCond = rgbMat(i,3) < 70;
    rgbCond =rCond && gCond && bCond; 
    if (cMat(i,2) < 0.45) && (cMat(i,3) > 0.1 && cMat(i,3) < 0.31) && rgbCond
        trashBin = [trashBin, i];
        msg = sprintf('%s \n Cluster #%d trashbin', msg, i);
%         noTree = false;
    end  
    %si esta entre aquests 3 nivells i a dalt, es arbre
    if(cMat(i, 1) < 0.5 && cMat(i, 2) < 0.25 && cMat(i, 3)<0.35 && mMat(i,1) <= mid)% && noTree
        msg = sprintf('%s \n Cluster #%d high tree prob. 1', msg, i);
        treeClustH = [treeClustH, i];
    end
    %Class light post
    if (cMat(i,1) > 0.4 && cMat(i,1) < 0.6) && (cMat(i,2) < 0.15) && (cMat(i,3) > 0.4 && cMat(i,3) < 0.6)
        msg = sprintf('%s \n Cluster #%d light post', msg, i);
        lightPost = [lightPost, i];
    end
end
%% Plot of classifications
%Exclou cluster que diferexen massa

if length(roadClustH) > 1
    roadClustDiff = sum(maxRoadSV(1) - svRatio(roadClustH(roadClustH ~= maxRoadSV(2))))/(length(roadClustH)-1);
    %si la diferencia entre clusters es massa alta, exclou el maxim i es queda amb els altres .
    if(roadClustDiff > 0.2)
        roadClustH(roadClustH == maxRoadSV(2)) = [];
    end
end

roadImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
treeImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
tbImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
farImg = zeros(size(hsvImg, 1), size(hsvImg, 2));

%Junta tots els clusters de carretera 
for i = roadClustH    
    roadImg = roadImg | clustImgs{i};
end
roadB = roadImg;
se = strel('disk', 1);
roadB = imopen(roadB, se);
roadB = imclose(roadB, se);
roadB(1:round(mid/2), :) = 0;
%exclou la part superior
roadImg(1:mid, :) = 0;
%junta els d'arbre
for i = treeClustH
    treeImg = treeImg | clustImgs{i};
end
%paperera
for i=trashBin
    tbImg= tbImg|clustImgs{i};
end
%light post
for i=lightPost
    farImg= farImg|clustImgs{i};
end

%creates a flat disk-shaped structuring element with the specified radius,
%R. Per obrir i tancar
se = strel('disk', 1);
threshRoad = imopen(roadImg, se);
threshRoad = imclose(threshRoad, se);

se = strel('rectangle', [3 1]);
threshTree = imopen(treeImg, se);
se = strel('rectangle', [4 1]);
threshTree = imclose(threshTree, se);

se = strel('rectangle', [20 1]);
threshLP = imerode(farImg, se);
threshLP = imdilate(threshLP, se);

% se = strel('square', 3);
% threshTB = imerode(tbImg, se);
% tbImg = imdilate(threshTB, se);

%% Connected component & Polygon (Road)
if road
    CC = bwconncomp(threshRoad,4);
    roadList = CC.PixelIdxList;

    roadComponents = [];

    for i = 1:length(roadList)
        if(size(roadList{i}, 1) > 150)
            roadComponents = [roadComponents, roadList(i)];
        end
    end

    roadPolys = [];

    for i = 1:length(roadComponents)
        roadComp = false(size(hsvImg, 1), size(hsvImg, 2));
        roadComp(roadComponents{i}) = 1;

        route = FindEdge(roadComp);

        polyEdge = [route(:, 2), route(:, 1)];% o route(:,2:-1:1) es el mateix
        
        newEdge = PolygonReduction(polyEdge, 5, 0);

        roadPolys = [roadPolys; mat2cell(newEdge, size(newEdge, 1), size(newEdge, 2))];
    end
end
%% Connected component & Polygon (Tree)
if tree
    CC2 = bwconncomp(threshTree,4);
    [treeList, match] = checkSplitTrunk(CC2, 3, mid);
    treeComponents = [];
    %selecciona els que tenen prou pixels i forma d'arbre
    for i=1:length(treeList)
        if(size(treeList{i}, 1) > 30)
        [Y,X]=ind2sub(CC2.ImageSize,treeList{i});
        yEnds=[min(Y) max(Y)];
        xEnds=[min(X) max(X)];
        yEndsDiff=yEnds(2)-yEnds(1);
        xEndsDiff=xEnds(2)-xEnds(1);
            if mean(Y) <= mid %Descarta els que queden sota de la meitat de la imatge
                if size(treeList{i}, 1) > 100 && yEndsDiff/xEndsDiff >= 2
                    treeComponents = [treeComponents, treeList(i)];
                elseif size(treeList{i}, 1) <=100 && yEndsDiff/xEndsDiff > 3
                    treeComponents = [treeComponents, treeList(i)];
                end
            end
        end
    end

    treePolys = [];
    treeArea=[];
    cornersRectangle=[];
    centersRectangle=[];
    for i = 1:length(treeComponents)
        %for tree rectangle
        [Y,X] = ind2sub(size(hsvImg),treeComponents{i});   
        corners=ApproxRectangle([X,Y]);
        cornersRectangle=[cornersRectangle; mat2cell(corners, size(corners, 1), size(corners, 2))];
            %for storage, probably later check with odometry.
    %         centers=[(corners(1,1)+corners(3,1))/2 (corners(1,2)+corners(2,2))/2];
    %         centersRectangle=[centersRectangle ; centers];
    end
end
%% Check for red Sign
if redSign
    [redImg, threshPoleImg, signCoords, rsBox, minCol, poleLines, lineImg] = FindRed(hImg,hsvImg2,'sign');
end
% check for red mid sign
if midSign
    [rmImg, ~, ~, msBox, ~, lines, postBelow] = FindRed(hImg,hsvImg2,'mid');
end
if bigSign
    [rbImg, ~, ~, bsBox, ~, ~, ~] = FindRed(hImg,hsvImg2,'big');
end
%% Check for trashbin
if trashbin
    if ~isempty(trashBin)
       [threshTB,tbBox]=FindTrashBin(tbImg); 
    end
end
%% LightPost
if lightpost
    if ~isempty(lightPost)
       [lpBox]=FindLightPost(threshLP); 
    end
end
toc
%% Storage
% if first == true 
%     previousCenters = [0 0];
%     previousCorners = zeros(4,2);
%     previousAreas = 0;
%     first = true;
% else
%      previousCenters = currentCenters;
%      previousCorners = currentCorners;
%      previousAreas = currentAreas;
% end
% currentCenters = centersRectangle;
% currentCorners = cornersRectangle;
% currentAreas = treeArea;


%% Figures
disp(msg)
if plots == true
figure(1)
imshow(Img, [])
title('Original image');

figure(2)
imshow(sImg, [])
title('Scaled image');

splot=ceil(cluster/4);
subPlotSize = [splot, 4];
figure(3)
for cClust = 1:cluster
    subplot(subPlotSize(1, 1), subPlotSize(1, 2), cClust)
    imshow(clustImgs{cClust}, [])
    title(['Cluster: ', num2str(cClust)])
end

if road
figure(4)
imshow(threshRoad, [])
title('Road Image')

if ~isempty(roadComponents)
    fig5 = figure(5);
    imshow(sImg,[]);
%     handl2 = fig5.CurrentAxes;
    for i = 1:length(roadPolys)
        polygon = roadPolys{i};
        k = impoly(gca, polygon);
    end
end
end

if tree
        if ~isempty(treeClustH)
        figure(6)
        imshow(treeImg, [])
        title('Tree Image')
        end


    if ~isempty(treeClustH)
        figure(7)
        imshow(threshTree, [])
        title('Tree Processed')
    end
    
    if ~isempty(treeComponents)
    fig8 = figure(8);
    imshow(sImg,[]);
%     handl4 = fig8.CurrentAxes;
    for i = 1:length(cornersRectangle)
        k = impoly(gca, cornersRectangle{i});
    end
    end
end

if trashbin
if ~isempty(trashBin)
if ~isempty(tbBox)
   figure(9)
   imshow(sImg)
   imrect(gca,tbBox);
end
end
end

if redSign
if ~(isempty(threshPoleImg))
    if ~(isempty(redImg) || isempty(rsBox))
        figure(10)
        rImg=hImg;
        rImg(redImg)=255;
        imshow(rImg)
        imrect(gca,rsBox);
    end
    figure(11)
    imshow(threshPoleImg, [])
    

    y = 1:size(lineImg, 1);
    xL = poleLines(2,1)*y + poleLines(1,1);
    xR = poleLines(2,2)*y + poleLines(1,2);
    figure(12)
    imshow(lineImg, [])
    hold on
    plot(xL, y, 'LineWidth', 2)
    plot(xR, y, 'LineWidth', 2)
    hold off
end
end

if midSign
if ~isempty(msBox)
   figure(13)
    rImg=hImg;
    rImg(rmImg)=255;
    imshow(rImg)
    imrect(gca,msBox);
if isfield(lines, 'point1')
        figure(14), imshow(postBelow), hold on
        max_len = 0;
        for k = 1:length(lines)
           xy = [lines(k).point1; lines(k).point2];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

           % Plot beginnings and ends of lines
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

           % Determine the endpoints of the longest line segment
           len = norm(lines(k).point1 - lines(k).point2);
           if ( len > max_len)
              max_len = len;
              xy_long = xy;
           end
        end
end
end
end

if bigSign
    if ~isempty(bsBox)
    figure(15)
    rImg=hImg;
    rImg(rbImg)=255;
    imshow(rImg)
    imrect(gca,bsBox);
    end
end
if lightpost
if ~isempty(lightPost)
    if ~isempty(lpBox)
        figure(16)
        imshow(sImg)
        imrect(gca,lpBox);
    end
end
end

end