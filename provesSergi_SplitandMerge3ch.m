% clc, close all, clear all
clear all
smoothing=1;
plotting=0;
plotstep = 0;
thvec=[0.1,0.05,0.05];  % threshold for red,green and blue (values b/w 0 and 1)
picname='r2.bmp';
I = imread(picname);  
I = imresize(I, [256 256]);
RGB = I;
I = rgb2hsv(I); 
if smoothing
    I2=I;
    I = GaussSmooth(I,16,2);
end

minBlockSize = 1;
tic
classtr=class(I);
[ir ic id]=size(I); % size of image before padding
if(id~=3)
    error('Input Matrix (Image) must be of size M-by-N-by-3 (true color)');    
end
if(~ispowerof2(ir) || ~ispowerof2(ic))
    I=padrgbtomakepowof2(I);
end
[irp,icp,id]=size(RGB); % size of image after padding
s=zeros(irp,icp);
H = I(:,:,1);
S = I(:,:,2);
V = I(:,:,3);
% % Matrix S is intitially of same size as input image/matrix and holds
% % zeros in all locations.
% % stru1 holds information about those images that are to be tested
% % for splitting.
% % initially we would put coordinates of ULC of image into stru1 then
% % we will execute a loop until stru1 becomes empty. 
% % Inside loop we would take coordinates of last image block from stru1
% % and remove them from stru1, then test if for splitting.
% % If it split then we would put four new blocks into stru1. If it does
% % not split then we would (ULC) of that single block into matrix S.
% % This process would continue until stru1 becomes emtpy, that means no
% % more splitting is required and we would be outside while loop. Then
% % simple we would converte matrix S into sparse matrix.

count1=1;                    
stru1(count1).rc=[1 1];      
stru1(count1).sz=size(I,1);  
while(size(stru1,2)>=1)    
    r1=stru1(end).rc(1);
    r2=r1+stru1(end).sz-1;
    c1=stru1(end).rc(2);
    c2=c1+stru1(end).sz-1; 
    stru1(end)=[];          % remove from stru1
    count1=count1-1;
    % % splitting (return splitting coordiantes with respect to segment
    % % (block) of image (matrix)
    [b1r,b1c,b2r,b2c,b3r,b3c,b4r,b4c,sz]=qtrgbsplitSergi(I(r1:r2,c1:c2,1:3),minBlockSize,thvec);
    if(sz==-1) % Not splitted
        % % disp('not splitted');
        s(r1,c1)=r2-r1+1;
    else      % splitted 
        % % disp('splitted');
        % % converting local coordinates to global coordinates
        b1r=r1+b1r-1;
        b1c=c1+b1c-1;
        b2r=r1+b2r-1;
        b2c=c1+b2c-1;
        b3r=r1+b3r-1;
        b3c=c1+b3c-1;
        b4r=r1+b4r-1;
        b4c=c1+b4c-1;        
        % % adding fourth block data in stru1
        count1=count1+1;
        stru1(count1).rc=[b4r,b4c];      
        stru1(count1).sz=sz;
        % % adding third block data in stru1
        count1=count1+1;
        stru1(count1).rc=[b3r,b3c];      
        stru1(count1).sz=sz;             
        % % adding second block data in stru1
        count1=count1+1;
        stru1(count1).rc=[b2r,b2c];      
        stru1(count1).sz=sz;             
        % % adding first block data in stru1
        count1=count1+1;
        stru1(count1).rc=[b1r,b1c];      
        stru1(count1).sz=sz;
    end    
end
s=sparse(s); %converting S into sparse matrix
if plotting
figure,imshow(RGB), hold on
end
% LABELING & MEANS
[y,x,blkSz] = find(s);
blkLgth = blkSz .^2;
labels = zeros([size(I,1) size(I,2)]);
Hmeans = zeros(length(y),1);
Vmeans = zeros(length(y),1);
Smeans = zeros(length(y),1);
for label=1:length(y) 
        labels(y(label):y(label)+blkSz(label)-1,x(label):x(label)+blkSz(label)-1) = label;
        Hmeans(label)=mean2(H(y(label):y(label)+blkSz(label)-1,x(label):x(label)+blkSz(label)-1));
        Vmeans(label)=mean2(S(y(label):y(label)+blkSz(label)-1,x(label):x(label)+blkSz(label)-1));
        Smeans(label)=mean2(V(y(label):y(label)+blkSz(label)-1,x(label):x(label)+blkSz(label)-1));
        stru1(label).lb = label;
        if plotting
            xp = [x(label) x(label)+blkSz(label) x(label)+blkSz(label) x(label)];
            yp = [y(label) y(label) y(label)+blkSz(label) y(label)+blkSz(label)];
            plot(xp,yp,'-','Color','r')
            hold on
        end
end
toc
%%
tic
previousLength = 0;
currentLength = 1;
iteration = 0;
hThr = 0.1; sThr = 0.1; vThr = 0.05; 
% MERGE
while 1
[nodes neighbours] = imRAG(labels);
currentLength = length(neighbours);
if abs(previousLength-currentLength) < 100
    break;
end
previousLength = currentLength;
label1 = neighbours(:,1);
label2 = neighbours(:,2);
for i=1:length(label1)

    hCond = abs(Hmeans(label1(i)) - Hmeans(label2(i)));
    sCond = abs(Smeans(label1(i)) - Smeans(label2(i)));
    vCond = abs(Vmeans(label1(i)) - Vmeans(label2(i)));
    if (hCond <= hThr) && (sCond <= sThr) && (vCond <= vThr)
       for nn = 1:length(stru1(label2(i)).lb)
           labels(y(stru1(label2(i)).lb(nn)):y(stru1(label2(i)).lb(nn))+blkSz(stru1(label2(i)).lb(nn))-1,x(stru1(label2(i)).lb(nn)):x(stru1(label2(i)).lb(nn))+blkSz(stru1(label2(i)).lb(nn))-1) = label1(i);
       end
       label1(label1 == label2(i)) = label1(i);
       
       stru1(label1(i)).lb = [stru1(label1(i)).lb stru1(label2(i)).lb];
       stru1(label2(i)).lb = 0;
       
       ttlSz = blkLgth(label1(i)) + blkLgth(label2(i)); 
       Hmeans(label1(i)) = Hmeans(label1(i))*blkLgth(label1(i))/ttlSz + Hmeans(label2(i))*blkLgth(label2(i))/ttlSz;
       Smeans(label1(i)) = Smeans(label1(i))*blkLgth(label1(i))/ttlSz + Smeans(label2(i))*blkLgth(label2(i))/ttlSz;
       Vmeans(label1(i)) = Vmeans(label1(i))*blkLgth(label1(i))/ttlSz + Vmeans(label2(i))*blkLgth(label2(i))/ttlSz;
       Hmeans(label2(i)) = -1;
       Smeans(label2(i)) = -1;
       Vmeans(label2(i)) = -1;
       blkLgth(label1(i)) = ttlSz;
       blkLgth(label2(i)) = 0;
    end
end
%PLOT EVERY STEP
if plotstep
    minSizePlot = size(I,1)/2-1;
    figure,imshow(RGB)
    hold on;
    zeroim = zeros(size(I));
    for i = 1:max(label1)
        if ~isempty(find(labels == i, 1)) && length(find(labels == i)) > minSizePlot
            zeroim(labels == i) = 1;
            zeroim(labels ~= i) = 0;
            route = FindEdge(zeroim);
            polyEdge = [route(:, 2), route(:, 1)];
            newEdge = PolygonReduction(polyEdge, 1, 0);
            plot(newEdge(:,1), newEdge(:,2), 'linewidth', 1, 'color', 'r');
        end
    end
end
iteration = iteration + 1;
end
toc
%%
I=I(1:ir,1:ic,:);
labels=labels(1:ir,1:ic,:);
minSizePlot = round(size(I,1)*2);
cstring='rgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmyk';
        figure,imshow(RGB)
        hold on;
        color=0;
        numClust = 0;
zeroim = zeros(size(I,1),size(I,2));
for i = 1:max(label1)
    if ~isempty(find(labels == i, 1)) && length(find(labels == i)) > minSizePlot
        color=color+1;
        numClust =numClust +1;
        zeroim(labels == i) = 1;
        zeroim(labels ~= i) = 0;
        zeroim=bwareaopen(zeroim,4,4);
        route = FindEdge(zeroim);
        polyEdge = [route(:, 2), route(:, 1)];
        newEdge = PolygonReduction(polyEdge, 1, 0);
        plot(newEdge(:,1), newEdge(:,2), 'linewidth', 2, 'color', cstring(color));
    end
end

% set(gca,'LooseInset',get(gca,'TightInset'));
%% Plot of X label
plabel = 0;
histo = 0;
sample = 580;%label2(find(label1 == 10));
if plabel == true
    bw = zeros(size(labels));
    for i=1:length(sample)
        bw(labels==sample(i)) = 1;
    end
    figure, imshow(bw)
    title(['Mean: ', num2str(mean2(H((labels==sample)))),' ',num2str(mean2(S((labels==sample)))),' ',num2str(mean2(V((labels==sample))))])
end

if histo
    for i=1:length(sample)
       figure
       hist(H(labels == sample(i)),0.05:0.05:0.95)
       figure
       hist(S(labels == sample(i)),0.05:0.05:0.95)
       figure
       hist(V(labels == sample(i)),0.05:0.05:0.95)
    end
end    