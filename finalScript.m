clc
close all
clear all

%% Import images

tImg = imread('nc7.bmp');

Img = tImg(:, 1:end-4, :);

figure(1)
imshow(Img, [])
title('Original image');
size(Img);
scale = 1/5;
sImg = imresize(Img, scale);
size(sImg);
hImg = imresize(Img, 1/2);

figure(2)
imshow(sImg, [])
title('Scaled image');

%% HSV Tranformation
tic
hsvImg = rgb2hsv(sImg);
hsvImg2 = rgb2hsv(hImg);

%cada H S V els ordena per columnes
lenLayer = size(hsvImg,1)*size(hsvImg,2);
rawInput = [hsvImg(1:lenLayer)' hsvImg((lenLayer + 1):(lenLayer*2))' hsvImg((lenLayer*2 + 1):(lenLayer*3))'];

%% 4D Histogram

%Input és la Imatge i l'interval
[hist3D, space, pixCor] = Hist3D(hsvImg, 25);
%
%% Seeding

clustThresh = 50;
clustCoords = [];
clustVal = clustThresh + 1;
oldClustVal = 0;
count = 1;

peaks3D = FindProminentPeak4D(hist3D, 1, 2);

while(true)
    
    if(count == 1)
        %com ones(), pero el valor és lògic al workspace
        peakMask = true(length(peaks3D) , 1);
        
    else
        peakMask = peaks3D(:, 1) > (coords(1, 1) + 1) | peaks3D(:, 1) < (coords(1, 1) - 1) & ...
                   peaks3D(:, 2) > (coords(1, 2) + 1) | peaks3D(:, 2) < (coords(1, 2) - 1) & ...
                   peaks3D(:, 3) > (coords(1, 3) + 1) | peaks3D(:, 3) < (coords(1, 3) - 1);
               
%         peakMask = (peaks3D(:, 1) > (coords(1, 1) + 1) | peaks3D(:, 1) < (coords(1, 1) - 1) )| ...
%                    (peaks3D(:, 2) > (coords(1, 2) + 1) | peaks3D(:, 2) < (coords(1, 2) - 1) )| ...
%                    (peaks3D(:, 3) > (coords(1, 3) + 1) | peaks3D(:, 3) < (coords(1, 3) - 1));

        peakMask = peakMask & oldMask;
    end
    
    newPeaks = peaks3D(peakMask, :);
    
    clustVal = max(newPeaks(:,4));

    if(clustVal < clustThresh)
        break;
    end
    %Enganxa les coordenades del colega amb valor més alt
    coords = newPeaks(newPeaks(:,4) == max(newPeaks(:,4)), 1:3);
    
    clustCoords = [clustCoords; coords];
    
    oldMask = peakMask;
    
    oldClustVal = clustVal;
    
    count = count + 1;
end

%% K-Means clustering

seed = clustCoords;
cluster = size(seed, 1);

[idx, cMat, sumd] = kmeans(rawInput, cluster, 'Start', seed./25);

%% Plotting of clusters

[hx, hy] = ind2sub([size(hsvImg, 1) size(hsvImg, 2)], 1:lenLayer);

figure(5)

mMat = zeros(cluster, 2);

splot=ceil(cluster/4);
% subPlotSize = factor(cluster + 1);
% if(length(subPlotSize) == 1)
%     subPlotSize = factor(cluster + 1 + 1);
% end
subPlotSize = [splot, 4];

clustImgs = cell(1, cluster);

for cClust = 1:cluster
    clustImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
    %Extreu X i Y de cada punt/píxel del cluster
    cXVal = hx(idx == cClust);
    cYVal = hy(idx == cClust);
    %assigna 1 al cluster 0 a la resta, per la binarització de la gràfica
    for i = 1:length(cXVal)
        clustImg(cXVal(i), cYVal(i)) = 1;
    end
    %Guarda el clúster a la seva corresponent posició de la cel·la    
    clustImgs{cClust} = clustImg;
    %Coordenada mitja dels pixels en la imatge
    mMat(cClust, :) = [round(sum(cXVal)/length(cXVal)), round(sum(cYVal)/length(cYVal))];
    
    subplot(subPlotSize(1, 1), subPlotSize(1, 2), cClust)
    imshow(clustImg, [])
    title(['Cluster: ', num2str(cClust)])
end
% subplot(subPlotSize(1, 1), subPlotSize(1, 2), cClust+1)
% imshow(sImg, [])
% title('Original Image');

%% Classification
 
clusters = 1:cluster;

cMat; % H S i V de cada centre de clúster

mid = size(hsvImg, 1)/2;
%Associa el cel al V més alt de tots
maxSkyClust = find(cMat(:, 3) == max(cMat(:, 3)));

skyClust = [];
%Si la saturació esta per sota de 0.1 i la mitjana de pixels per sobre la
%meitat de la foto, l'assigna i l'esborra.
if(cMat(maxSkyClust, 2) < 0.1 && mMat(maxSkyClust, 1) < mid)
    skyClust = maxSkyClust;
    strcat('The sky cluster is: #', num2str(maxSkyClust))
    clusters(skyClust) = [];
end
%Calcula ratios
hvRatio = cMat(:, 1) ./ cMat(:, 3);
svRatio = cMat(:, 2) ./ cMat(:, 3);
hsRatio = cMat(:, 1) ./ cMat(:, 2);

% [hsRatio, hvRatio, svRatio];

roadClustH = [];
roadClustL = [];
treeClustH = [];
treeClustL = [];
maxRoadSV = zeros(1,2);
msg = '';
for i = clusters
    foundTreeClust = 0;
    %usa un ratio i un màx de value per la carretera
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
            %els que no han passat el tall pero compleixen aquestes
            %condicions seran carretera de poca probabilitat.
        elseif(cMat(i, 3) > 0.30 && mMat(i, 1) >= mid - 15)
            roadClustL = [roadClustL, i];
            msg = sprintf('%s \n Cluster #%d low\t road prob.', msg, i);
        end
    end
    %si esta entre aquests 3 nivells i no sha trobat abans, es arbre
    if(cMat(i, 1) < 0.26 && cMat(i, 1) > 0.17 && ...
       cMat(i, 2) < 0.24 && cMat(i, 2) > 0.16 && ...
       cMat(i, 3) < 0.37 && cMat(i, 3) > 0.28 && ...
       foundTreeClust == 0)
        msg = sprintf('%s \n Cluster #%d high tree prob. 1', msg, i);
        treeClustH = [treeClustH, i];
        foundTreeClust = 1;
    end
    %no entence la condicio de foundtree clust
    if(cMat(i, 1) < 0.30 && cMat(i, 1) > 0.13 && ...
       cMat(i, 2) < 0.28 && cMat(i, 2) > 0.12 && ...
       cMat(i, 3) < 0.41 && cMat(i, 3) > 0.24 && ...
       foundTreeClust == 0)
        msg = sprintf('%s \n Cluster #%d low\t tree prob.', msg, i);
        treeClustL = [treeClustL, i];
        foundTreeClust = 1;
    end
    
    %     if(hvRatio(i) < 0.85 && hvRatio(i) > 0.60 && ...
%        svRatio(i) < 0.70 && svRatio(i) > 0.50 && ...
%        hsRatio(i) > 0.70 && foundTreeClust == 0)
%         msg = sprintf('%s \n Cluster #%d high tree prob. 2', msg, i);
%         treeClustH = [treeClustH, i];
%         foundTreeClust = 1;
%     end
end

msg

roadImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
treeImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
treeImgL = zeros(size(hsvImg, 1), size(hsvImg, 2));

roadClustDiff = sum(maxRoadSV(1) - svRatio(roadClustH(roadClustH ~= maxRoadSV(2))))/(length(roadClustH)-1);
%si la diferencia entre clusters es massa alta, exclou els altres i es
%queda amb el maxim.
if(roadClustDiff > 0.2)
    roadClustH(roadClustH == maxRoadSV(2)) = [];
end
%Junta tots els clusters de carretera 
for i = roadClustH    
    roadImg = roadImg | clustImgs{i};
end
%exclou la part superior
roadImg(1:mid, :) = 0;
%junta els d'arbre
for i = treeClustH
    treeImg = treeImg | clustImgs{i};
end
for i = treeClustL
    treeImgL = treeImgL | clustImgs{i};
end
%creates a flat disk-shaped structuring element with the specified radius,
%R. Per obrir i tancar
se = strel('disk', 1);
threshRoad = imopen(roadImg, se);
threshRoad2 = imclose(threshRoad, se);

figure(7)
imshow(threshRoad2, [])
title('Road Image')
%%
if ~isempty(treeClustH)
figure(8)
imshow(treeImg, [])
title('Tree Image')
end

if ~isempty(treeClustL)
figure(11)
imshow(treeImgL, [])
title('Low probability Tree Image')
end
%% Connected component & Polygon (Road)

CC = bwconncomp(threshRoad2);
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
%     size(polyEdge)
    newEdge = PolygonReduction(polyEdge, 5, 0);

    roadPolys = [roadPolys; mat2cell(newEdge, size(newEdge, 1), size(newEdge, 2))];
end
toc

if ~isempty(roadComponents)
    routeImg = false(size(roadComp));

%     fig8 = figure(9);
%     imshow(routeImg, []);
    fig12 = figure(12);
    imshow(sImg,[]);
%     handl = fig8.CurrentAxes;
    handl2 = fig12.CurrentAxes;
    for i = 1:length(roadPolys)
%         h = impoly(handl, roadPolys{i});
        k = impoly(handl2, roadPolys{i});
    end
end
% 
% %% Connected component (Trees)
% 
% CC = bwconncomp(treeImg);
% treeList = CC.PixelIdxList;
% 
% treeComponents = [];
% 
% threshTreeImg = zeros(size(hsvImg, 1), size(hsvImg, 2));
% for i = 1:length(treeList)
%     
%     if(size(treeList{i}, 1) > 200)
%         treeComponents = [treeComponents, treeList(i)];
%         threshTreeImg(treeList{i}) = 1;
%     end
% end
% 
% figure(10)
% imshow(threshTreeImg, [])

%% Check for landmark

[threshPoleImg, signCoords] = FindRedSign(hsvImg2);

if(~isempty(threshPoleImg))
    poleLines = FindParallelLines(threshPoleImg, signCoords);
else
    'There is no sign in the image'
end